// CLI entry point. cli/ is the top of the dependency graph (CLAUDE.md #3):
// it is the only place allowed to link core/, forcefield/, engines/, and
// io/ together.
//
// Design goals for this milestone (the CLI is meant to be a real UX
// differentiator, not just a wiring proof): config validation happens
// entirely up front — a malformed config fails immediately with a message
// naming the offending key (and, where toml++ can supply one, the exact
// line/column), never partway through a run. `--dry-run` runs every check
// a real run would (including opening and parsing the structure file) and
// prints the fully-resolved configuration without touching the physics.
// `--json` switches the final result to a single line of machine-readable
// JSON on stdout; progress and diagnostic messages always go to stderr via
// spdlog, so piping `--json` output never mixes the two.
//
// Exit codes: 0 = success. 1 = configuration/validation error (including
// an unrecognized command). 2 = the requested physics is genuinely not
// implemented yet (engines/dynamics, engines/geometry_analysis) — surfaced
// as a clean, expected NotImplemented rather than a crash or silent no-op,
// but distinguished from a config error so a calling script can tell the
// two apart.

#include <algorithm>
#include <array>
#include <chrono>
#include <cstdlib>
#include <exception>
#include <filesystem>
#include <iostream>
#include <memory>
#include <sstream>
#include <string>
#include <string_view>
#include <vector>

#include <aleator/version.hpp>
#include <spdlog/spdlog.h>

#include "core/exceptions.hpp"
#include "core/geometry/lattice.hpp"
#include "core/math/counter_based_rng.hpp"
#include "core/math/particle_data.hpp"
#include "core/neighbor/verlet_list.hpp"
#include "engines/dynamics/integrator.hpp"
#include "engines/geometry_analysis/pore_analysis.hpp"
#include "engines/monte_carlo/molecule_species.hpp"
#include "engines/monte_carlo/monte_carlo_engine.hpp"
#include "engines/monte_carlo/peng_robinson.hpp"
#include "forcefield/pairwise/lennard_jones.hpp"
#include "io/config.hpp"
#include "io/config_error.hpp"
#include "io/structure_io.hpp"

namespace {

constexpr int kExitSuccess = EXIT_SUCCESS;
constexpr int kExitConfigError = 1;
constexpr int kExitNotImplemented = 2;

void printUsage() {
    std::cout
        << "aleator " << aleator::kVersion << "\n\n"
        << "Usage:\n"
        << "  aleator --version                              Print the version and exit\n"
        << "  aleator --help                                 Print this message and exit\n"
        << "  aleator gcmc run <config.toml> [flags]          Run grand-canonical Monte Carlo\n"
        << "  aleator pore analyze <config.toml> [flags]      Run pore-geometry analysis\n"
        << "  aleator md run <config.toml> [flags]            Run molecular dynamics\n"
        << "  aleator validate <config.toml>                  Validate a config, then exit\n"
        << "  aleator bench [--json]                          Run a quick built-in benchmark\n"
        << "\n"
        << "Flags (gcmc/pore/md):\n"
        << "  --dry-run   Validate the config and structure file, print the resolved\n"
        << "              configuration, and exit without running anything.\n"
        << "  --json      Emit the final result as one line of JSON on stdout instead of\n"
        << "              human-readable text. Progress/diagnostic messages still go to\n"
        << "              stderr, so this is safe to pipe.\n"
        << "\n"
        << "Exit codes: 0 success, 1 configuration/usage error, 2 requested feature not\n"
        << "implemented yet in this build.\n";
}

std::string jsonEscape(std::string_view s) {
    std::string out;
    out.reserve(s.size());
    for (char c : s) {
        switch (c) {
            case '"':
                out += "\\\"";
                break;
            case '\\':
                out += "\\\\";
                break;
            case '\n':
                out += "\\n";
                break;
            default:
                out += c;
        }
    }
    return out;
}

struct ParsedArgs {
    std::vector<std::string_view> positional;
    bool dryRun = false;
    bool json = false;
};

ParsedArgs parseFlags(const std::vector<std::string_view>& args) {
    ParsedArgs result;
    for (auto arg : args) {
        if (arg == "--dry-run") {
            result.dryRun = true;
        } else if (arg == "--json") {
            result.json = true;
        } else {
            result.positional.push_back(arg);
        }
    }
    return result;
}

// ---------------------------------------------------------------- gcmc ---

struct PreparedGcmc {
    aleator::io::StructureData structure;
    std::shared_ptr<const aleator::forcefield::LennardJones> forceField;
    aleator::engines::MoleculeSpecies adsorbate;
    double fugacityPascal;
    bool idealGas;
};

/// Loads the framework structure and builds the LJ force field + adsorbate
/// species template — everything a real run needs, and everything
/// --dry-run should validate but not act on. Throws if the CIF is
/// malformed or if `gcmc.framework_lj` doesn't cover every element the CIF
/// actually contains (a config that's silently missing a species is
/// exactly the kind of failure that should surface at startup, not as a
/// wrong energy mid-run).
PreparedGcmc prepareGcmc(const aleator::io::GcmcRunConfig& cfg) {
    auto structure = aleator::io::readCif(cfg.gcmc.frameworkCif);

    std::vector<aleator::forcefield::LennardJonesParameters> ljParameters(
        structure.speciesSymbols.size());
    for (std::size_t i = 0; i < structure.speciesSymbols.size(); ++i) {
        const std::string& symbol = structure.speciesSymbols[i];
        const auto it = std::find_if(cfg.gcmc.frameworkLennardJones.begin(),
                                      cfg.gcmc.frameworkLennardJones.end(),
                                      [&](const auto& s) { return s.element == symbol; });
        if (it == cfg.gcmc.frameworkLennardJones.end()) {
            throw std::runtime_error("gcmc.framework_lj has no entry for element \"" + symbol +
                                      "\", which is present in " +
                                      cfg.gcmc.frameworkCif.string() +
                                      " -- add a [[gcmc.framework_lj]] entry for it");
        }
        ljParameters[i] = {it->epsilonKelvin, it->sigmaAngstrom};
    }
    const auto adsorbateSpeciesIndex = static_cast<std::uint32_t>(ljParameters.size());
    ljParameters.push_back({cfg.gcmc.adsorbate.epsilonKelvin, cfg.gcmc.adsorbate.sigmaAngstrom});

    auto forceField = std::make_shared<const aleator::forcefield::LennardJones>(
        ljParameters, cfg.gcmc.cutoffAngstrom, aleator::forcefield::LennardJonesTruncation::Shifted);

    aleator::engines::MoleculeSpecies adsorbate;
    adsorbate.sites.push_back({0.0, 0.0, 0.0, adsorbateSpeciesIndex, cfg.gcmc.adsorbate.chargeE,
                                cfg.gcmc.adsorbate.massAmu});

    const double pressurePascal = cfg.gcmc.pressureBar * 1.0e5;
    double fugacityPascal = pressurePascal;
    bool idealGas = true;
    if (cfg.gcmc.adsorbate.criticalProperties.has_value()) {
        const auto& cp = *cfg.gcmc.adsorbate.criticalProperties;
        const aleator::engines::PengRobinson eos(
            {cp.criticalTemperatureKelvin, cp.criticalPressurePascal, cp.acentricFactor});
        fugacityPascal = eos.fugacityPascal(cfg.gcmc.temperatureKelvin, pressurePascal);
        idealGas = false;
    }

    return {std::move(structure), std::move(forceField), std::move(adsorbate), fugacityPascal,
            idealGas};
}

void printGcmcConfig(const aleator::io::GcmcRunConfig& cfg, bool json) {
    const auto& run = cfg.run;
    const auto& g = cfg.gcmc;
    if (json) {
        std::ostringstream oss;
        oss << "{\"run\":{\"name\":\"" << jsonEscape(run.name) << "\",\"rng_seed\":" << run.rngSeed
            << ",\"thread_count\":" << run.threadCount << "},"
            << "\"gcmc\":{\"framework_cif\":\"" << jsonEscape(g.frameworkCif.string())
            << "\",\"temperature_kelvin\":" << g.temperatureKelvin
            << ",\"pressure_bar\":" << g.pressureBar << ",\"cutoff_angstrom\":" << g.cutoffAngstrom
            << ",\"equilibration_steps\":" << g.equilibrationSteps
            << ",\"production_steps\":" << g.productionSteps << ",\"adsorbate\":{\"name\":\""
            << jsonEscape(g.adsorbate.name) << "\",\"epsilon_kelvin\":" << g.adsorbate.epsilonKelvin
            << ",\"sigma_angstrom\":" << g.adsorbate.sigmaAngstrom
            << ",\"mass_amu\":" << g.adsorbate.massAmu << ",\"charge_e\":" << g.adsorbate.chargeE
            << ",\"uses_peng_robinson\":"
            << (g.adsorbate.criticalProperties.has_value() ? "true" : "false")
            << "},\"framework_lj\":[";
        for (std::size_t i = 0; i < g.frameworkLennardJones.size(); ++i) {
            const auto& s = g.frameworkLennardJones[i];
            if (i > 0) {
                oss << ",";
            }
            oss << "{\"element\":\"" << jsonEscape(s.element)
                << "\",\"epsilon_kelvin\":" << s.epsilonKelvin
                << ",\"sigma_angstrom\":" << s.sigmaAngstrom << "}";
        }
        oss << "]}}";
        std::cout << oss.str() << "\n";
        return;
    }

    std::cout << "run.name                      = " << run.name << "\n"
              << "run.rng_seed                  = " << run.rngSeed << "\n"
              << "run.thread_count               = " << run.threadCount << "\n"
              << "gcmc.framework_cif             = " << g.frameworkCif.string() << "\n"
              << "gcmc.temperature_kelvin        = " << g.temperatureKelvin << "\n"
              << "gcmc.pressure_bar              = " << g.pressureBar << "\n"
              << "gcmc.cutoff_angstrom           = " << g.cutoffAngstrom << "\n"
              << "gcmc.equilibration_steps       = " << g.equilibrationSteps << "\n"
              << "gcmc.production_steps          = " << g.productionSteps << "\n"
              << "gcmc.adsorbate.name            = " << g.adsorbate.name << "\n"
              << "gcmc.adsorbate.epsilon_kelvin  = " << g.adsorbate.epsilonKelvin << "\n"
              << "gcmc.adsorbate.sigma_angstrom  = " << g.adsorbate.sigmaAngstrom << "\n"
              << "gcmc.adsorbate.mass_amu        = " << g.adsorbate.massAmu << "\n"
              << "gcmc.adsorbate.charge_e        = " << g.adsorbate.chargeE << "\n"
              << "gcmc.adsorbate.fugacity_model  = "
              << (g.adsorbate.criticalProperties.has_value()
                      ? "Peng-Robinson"
                      : "ideal gas (no critical properties given -- fugacity = pressure)")
              << "\n"
              << "gcmc.framework_lj (" << g.frameworkLennardJones.size() << " species):\n";
    for (const auto& s : g.frameworkLennardJones) {
        std::cout << "  " << s.element << ": epsilon=" << s.epsilonKelvin
                  << " K, sigma=" << s.sigmaAngstrom << " Ang\n";
    }
}

int runGcmcCommand(const std::filesystem::path& configPath, bool dryRun, bool json) {
    const aleator::io::GcmcRunConfig cfg = aleator::io::loadGcmcConfig(configPath);

    if (dryRun) {
        prepareGcmc(cfg); // validates the CIF and species coverage too
        printGcmcConfig(cfg, json);
        return kExitSuccess;
    }

    PreparedGcmc prepared = prepareGcmc(cfg);
    spdlog::info("loaded {} framework atoms ({} species) from {}",
                 prepared.structure.particles.size(), prepared.structure.speciesSymbols.size(),
                 cfg.gcmc.frameworkCif.string());
    if (prepared.idealGas) {
        spdlog::info(
            "no Peng-Robinson critical properties given for \"{}\" -- using ideal-gas fugacity "
            "(fugacity = pressure)",
            cfg.gcmc.adsorbate.name);
    } else {
        spdlog::info("fugacity({} K, {} bar) = {:.3f} Pa (Peng-Robinson)", cfg.gcmc.temperatureKelvin,
                     cfg.gcmc.pressureBar, prepared.fugacityPascal);
    }

    auto rng = std::make_unique<aleator::core::Philox4x32Rng>();
    rng->seed(cfg.run.rngSeed, /*streamIndex=*/0);
    const double volume = std::abs(prepared.structure.lattice.volume());
    aleator::engines::MonteCarloEngine engine(std::move(prepared.structure.particles),
                                               prepared.structure.lattice, prepared.forceField,
                                               std::move(rng), cfg.gcmc.temperatureKelvin,
                                               prepared.adsorbate, prepared.fugacityPascal);

    spdlog::info("equilibrating ({} steps)...", cfg.gcmc.equilibrationSteps);
    engine.run(cfg.gcmc.equilibrationSteps);

    spdlog::info("production ({} steps)...", cfg.gcmc.productionSteps);
    constexpr int kProgressUpdates = 10;
    const std::size_t stepsPerUpdate =
        std::max<std::size_t>(1, cfg.gcmc.productionSteps / kProgressUpdates);
    double sumN = 0.0;
    std::size_t samples = 0;
    std::size_t stepsRun = 0;
    while (stepsRun < cfg.gcmc.productionSteps) {
        const std::size_t thisBatch = std::min(stepsPerUpdate, cfg.gcmc.productionSteps - stepsRun);
        engine.run(thisBatch);
        stepsRun += thisBatch;
        sumN += static_cast<double>(engine.moleculeCount()) * static_cast<double>(thisBatch);
        samples += thisBatch;
        spdlog::info("  [{:3.0f}%] step {}/{}  <N> so far = {:.3f}  molecules = {}",
                     100.0 * static_cast<double>(stepsRun) /
                         static_cast<double>(cfg.gcmc.productionSteps),
                     stepsRun, cfg.gcmc.productionSteps,
                     sumN / static_cast<double>(samples), engine.moleculeCount());
    }
    const double meanN = sumN / static_cast<double>(samples);
    const double loadingMoleculesPerAngstrom3 = meanN / volume;

    if (json) {
        std::ostringstream oss;
        oss << "{\"mean_molecule_count\":" << meanN << ",\"loading_molecules_per_angstrom3\":"
            << loadingMoleculesPerAngstrom3 << ",\"cell_volume_angstrom3\":" << volume << "}";
        std::cout << oss.str() << "\n";
    } else {
        std::cout << "Final result:\n"
                  << "  <N> (mean molecule count over production) = " << meanN << "\n"
                  << "  loading                                   = " << loadingMoleculesPerAngstrom3
                  << " molecules/Ang^3\n";
    }
    return kExitSuccess;
}

// ---------------------------------------------------------------- pore ---

void printPoreConfig(const aleator::io::PoreRunConfig& cfg, bool json) {
    if (json) {
        std::ostringstream oss;
        oss << "{\"run\":{\"name\":\"" << jsonEscape(cfg.run.name) << "\"},"
            << "\"pore\":{\"framework_cif\":\"" << jsonEscape(cfg.pore.frameworkCif.string())
            << "\",\"probe_radius_angstrom\":" << cfg.pore.probeRadiusAngstrom << "}}";
        std::cout << oss.str() << "\n";
        return;
    }
    std::cout << "run.name                    = " << cfg.run.name << "\n"
              << "pore.framework_cif          = " << cfg.pore.frameworkCif.string() << "\n"
              << "pore.probe_radius_angstrom  = " << cfg.pore.probeRadiusAngstrom << "\n";
}

int runPoreCommand(const std::filesystem::path& configPath, bool dryRun, bool json) {
    const aleator::io::PoreRunConfig cfg = aleator::io::loadPoreConfig(configPath);
    const auto structure = aleator::io::readCif(cfg.pore.frameworkCif); // validates either way

    if (dryRun) {
        printPoreConfig(cfg, json);
        return kExitSuccess;
    }

    try {
        const auto geometry =
            aleator::engines::analyzePoreGeometry(structure.particles, structure.lattice);
        (void)geometry;
    } catch (const aleator::NotImplemented& err) {
        spdlog::warn("{} -- pore geometry analysis is not implemented in this build yet",
                     err.what());
        return kExitNotImplemented;
    }
    return kExitSuccess;
}

// ------------------------------------------------------------------ md ---

void printMdConfig(const aleator::io::MdRunConfig& cfg, bool json) {
    if (json) {
        std::ostringstream oss;
        oss << "{\"run\":{\"name\":\"" << jsonEscape(cfg.run.name) << "\"},"
            << "\"md\":{\"structure_file\":\"" << jsonEscape(cfg.md.structureFile.string())
            << "\",\"timestep_ps\":" << cfg.md.timestepPicoseconds
            << ",\"num_steps\":" << cfg.md.numSteps << "}}";
        std::cout << oss.str() << "\n";
        return;
    }
    std::cout << "run.name             = " << cfg.run.name << "\n"
              << "md.structure_file    = " << cfg.md.structureFile.string() << "\n"
              << "md.timestep_ps       = " << cfg.md.timestepPicoseconds << "\n"
              << "md.num_steps         = " << cfg.md.numSteps << "\n";
}

int runMdCommand(const std::filesystem::path& configPath, bool dryRun, bool json) {
    const aleator::io::MdRunConfig cfg = aleator::io::loadMdConfig(configPath);
    auto structure = aleator::io::readCif(cfg.md.structureFile); // validates either way

    if (dryRun) {
        printMdConfig(cfg, json);
        return kExitSuccess;
    }

    // No force-field selection in the MD config schema yet (future scope,
    // see io/config.hpp's MdConfig doc comment) -- an empty species table
    // is enough to prove the wiring, matching the pattern already
    // established for pore/. The integrator itself is what's actually
    // being exercised here, and it throws NotImplemented regardless.
    auto forceField = std::make_shared<const aleator::forcefield::LennardJones>(
        std::vector<aleator::forcefield::LennardJonesParameters>{}, 12.0);
    aleator::engines::VelocityVerletIntegrator integrator;
    try {
        integrator.step(structure.particles, structure.lattice, *forceField,
                         cfg.md.timestepPicoseconds);
    } catch (const aleator::NotImplemented& err) {
        spdlog::warn("{} -- MD integration is not implemented in this build yet", err.what());
        return kExitNotImplemented;
    }
    return kExitSuccess;
}

// ------------------------------------------------------------ validate ---

int runValidateCommand(const std::filesystem::path& configPath) {
    const auto kind = aleator::io::detectConfigKind(configPath);
    switch (kind) {
        case aleator::io::ConfigKind::Gcmc: {
            const auto cfg = aleator::io::loadGcmcConfig(configPath);
            const auto prepared = prepareGcmc(cfg);
            std::cout << configPath.string() << ": valid GCMC config (\"" << cfg.run.name << "\", "
                      << prepared.structure.particles.size() << " framework atoms, adsorbate=\""
                      << cfg.gcmc.adsorbate.name << "\")\n";
            return kExitSuccess;
        }
        case aleator::io::ConfigKind::Pore: {
            const auto cfg = aleator::io::loadPoreConfig(configPath);
            const auto structure = aleator::io::readCif(cfg.pore.frameworkCif);
            std::cout << configPath.string() << ": valid pore-analysis config (\"" << cfg.run.name
                      << "\", " << structure.particles.size() << " framework atoms)\n";
            return kExitSuccess;
        }
        case aleator::io::ConfigKind::Md: {
            const auto cfg = aleator::io::loadMdConfig(configPath);
            std::cout << configPath.string() << ": valid MD config (\"" << cfg.run.name << "\")\n";
            return kExitSuccess;
        }
        case aleator::io::ConfigKind::Unknown:
            std::cerr << "error: " << configPath.string()
                      << " must have exactly one of [gcmc], [pore], or [md] at the top level\n";
            return kExitConfigError;
    }
    return kExitConfigError;
}

// --------------------------------------------------------------- bench ---

int runBenchCommand(bool json) {
    constexpr std::size_t kParticleCount = 4000;
    // Cutoff kept well under box size (~1/10) so the cell list actually
    // gets enough cells per axis to prune candidate pairs -- a cutoff
    // close to L/2 degenerates every cell's neighborhood to the whole box
    // (no pruning at all), which is a benchmark-parameter mistake, not a
    // CellList correctness issue (confirmed: an earlier version of this
    // function used a 40 Ang box with a 10 Ang cutoff, giving only 3 cells
    // per axis -- meaning literally every cell was in the wraparound
    // neighborhood of every other cell -- and ran for minutes).
    constexpr double kBoxLength = 60.0;
    constexpr double kCutoff = 6.0;
    constexpr int kRepeats = 20;

    std::array<std::array<double, 3>, 3> matrix{};
    matrix[0] = {kBoxLength, 0.0, 0.0};
    matrix[1] = {0.0, kBoxLength, 0.0};
    matrix[2] = {0.0, 0.0, kBoxLength};
    const aleator::core::Lattice lattice(matrix);

    aleator::core::ParticleData particles;
    particles.resize(kParticleCount);
    aleator::core::Philox4x32Rng rng;
    rng.seed(1, 0);
    for (std::size_t i = 0; i < kParticleCount; ++i) {
        particles.x[i] = rng.nextUniform() * kBoxLength;
        particles.y[i] = rng.nextUniform() * kBoxLength;
        particles.z[i] = rng.nextUniform() * kBoxLength;
        particles.species[i] = 0;
    }

    const std::vector<aleator::forcefield::LennardJonesParameters> species{{100.0, 3.0}};
    const aleator::forcefield::LennardJones lj(species, kCutoff,
                                                aleator::forcefield::LennardJonesTruncation::Shifted);

    aleator::core::VerletList neighbors(/*skinDistance=*/0.5);
    const auto buildStart = std::chrono::steady_clock::now();
    neighbors.build(particles, lattice, kCutoff);
    const auto buildEnd = std::chrono::steady_clock::now();

    double energy = 0.0; // accumulated so the compiler can't optimize the loop away
    const auto energyStart = std::chrono::steady_clock::now();
    for (int i = 0; i < kRepeats; ++i) {
        energy += lj.computeEnergy(particles, lattice, neighbors);
    }
    const auto energyEnd = std::chrono::steady_clock::now();

    const double buildMs = std::chrono::duration<double, std::milli>(buildEnd - buildStart).count();
    const double energyMsPerCall =
        std::chrono::duration<double, std::milli>(energyEnd - energyStart).count() / kRepeats;
    const std::size_t pairCount = neighbors.pairs().size();

    if (json) {
        std::ostringstream oss;
        oss << "{\"particle_count\":" << kParticleCount << ",\"pair_count\":" << pairCount
            << ",\"neighbor_list_build_ms\":" << buildMs
            << ",\"lj_energy_ms_per_call\":" << energyMsPerCall
            << ",\"energy_reduction_check\":" << energy << "}";
        std::cout << oss.str() << "\n";
    } else {
        std::cout << "aleator bench -- quick built-in check (this build's flags, single-threaded)\n"
                  << "  particles                    = " << kParticleCount << "\n"
                  << "  neighbor pairs                = " << pairCount << "\n"
                  << "  VerletList::build              = " << buildMs << " ms\n"
                  << "  LennardJones::computeEnergy    = " << energyMsPerCall << " ms/call (avg over "
                  << kRepeats << " calls)\n\n"
                  << "For the full benchmark suite with baseline/regression tracking, use the\n"
                  << "`bench` CMake preset instead:\n"
                  << "  cmake --preset bench && cmake --build --preset bench && "
                     "./scripts/run_benchmarks.sh\n";
    }
    return kExitSuccess;
}

} // namespace

int main(int argc, char** argv) {
    const std::vector<std::string_view> rawArgs(argv + 1, argv + argc);

    if (rawArgs.empty() || rawArgs.front() == "--help" || rawArgs.front() == "-h") {
        printUsage();
        return rawArgs.empty() ? kExitConfigError : kExitSuccess;
    }
    if (rawArgs.front() == "--version") {
        std::cout << aleator::kVersion << "\n";
        return kExitSuccess;
    }

    try {
        if (rawArgs.front() == "gcmc") {
            if (rawArgs.size() < 2 || rawArgs[1] != "run") {
                std::cerr << "error: usage: aleator gcmc run <config.toml> [--dry-run] [--json]\n";
                return kExitConfigError;
            }
            const auto parsed =
                parseFlags(std::vector<std::string_view>(rawArgs.begin() + 2, rawArgs.end()));
            if (parsed.positional.size() != 1) {
                std::cerr << "error: usage: aleator gcmc run <config.toml> [--dry-run] [--json]\n";
                return kExitConfigError;
            }
            return runGcmcCommand(std::string(parsed.positional[0]), parsed.dryRun, parsed.json);
        }
        if (rawArgs.front() == "pore") {
            if (rawArgs.size() < 2 || rawArgs[1] != "analyze") {
                std::cerr
                    << "error: usage: aleator pore analyze <config.toml> [--dry-run] [--json]\n";
                return kExitConfigError;
            }
            const auto parsed =
                parseFlags(std::vector<std::string_view>(rawArgs.begin() + 2, rawArgs.end()));
            if (parsed.positional.size() != 1) {
                std::cerr
                    << "error: usage: aleator pore analyze <config.toml> [--dry-run] [--json]\n";
                return kExitConfigError;
            }
            return runPoreCommand(std::string(parsed.positional[0]), parsed.dryRun, parsed.json);
        }
        if (rawArgs.front() == "md") {
            if (rawArgs.size() < 2 || rawArgs[1] != "run") {
                std::cerr << "error: usage: aleator md run <config.toml> [--dry-run] [--json]\n";
                return kExitConfigError;
            }
            const auto parsed =
                parseFlags(std::vector<std::string_view>(rawArgs.begin() + 2, rawArgs.end()));
            if (parsed.positional.size() != 1) {
                std::cerr << "error: usage: aleator md run <config.toml> [--dry-run] [--json]\n";
                return kExitConfigError;
            }
            return runMdCommand(std::string(parsed.positional[0]), parsed.dryRun, parsed.json);
        }
        if (rawArgs.front() == "validate") {
            if (rawArgs.size() != 2) {
                std::cerr << "error: usage: aleator validate <config.toml>\n";
                return kExitConfigError;
            }
            return runValidateCommand(std::string(rawArgs[1]));
        }
        if (rawArgs.front() == "bench") {
            const auto parsed =
                parseFlags(std::vector<std::string_view>(rawArgs.begin() + 1, rawArgs.end()));
            return runBenchCommand(parsed.json);
        }
    } catch (const aleator::io::ConfigError& err) {
        std::cerr << "error: " << err.what() << "\n";
        return kExitConfigError;
    } catch (const std::exception& err) {
        std::cerr << "error: " << err.what() << "\n";
        return kExitConfigError;
    }

    std::cerr << "error: unrecognized command \"" << rawArgs.front() << "\"\n";
    printUsage();
    return kExitConfigError;
}
