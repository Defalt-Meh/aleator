// CLI entry point. cli/ is the top of the dependency graph (CLAUDE.md #3):
// it is the only place allowed to link core/, forcefield/, engines/, and
// io/ together.
//
// CLAUDE.md milestone: "the CLI is the differentiator... expose ONLY
// subcommands backed by working engines." `gcmc run`, `pore analyze`,
// `validate`, and `bench` exist as executable subcommands — `md run` does
// not, because engines/dynamics is NotImplemented, and a subcommand that
// prints NotImplemented when run advertises capability this build doesn't
// have (worse than omitting it). The `[md]` config *schema* still exists
// (io/config.hpp) and `aleator validate` still recognizes it — schema
// validation is real and honest regardless of whether an engine exists yet
// — but it says plainly that the engine isn't implemented rather than
// silently calling the config "valid" in a way that could be mistaken for
// "runnable".
//
// Config validation happens entirely up front — a malformed config fails
// immediately with a message naming the offending key (and, where toml++
// can supply one, the exact line/column), never partway through a run.
// `--dry-run` runs every check a real run would (including opening and
// parsing the structure file) and prints the fully-resolved configuration,
// every default included, without touching the physics. `--json` switches
// the final result to a single line of machine-readable JSON on stdout;
// progress and diagnostic messages always go to stderr via spdlog, so
// piping `--json` output never mixes the two. A real (non-dry-run) `gcmc
// run` also writes the fully-resolved configuration (including the RNG
// seed) to `run.output_directory`, so a published result can be
// reproduced from its artifacts alone.
//
// Exit codes: 0 = success. 1 = configuration/validation error (including
// an unrecognized command). 2 = `aleator validate` was pointed at a
// syntactically valid `[pore]`/`[md]` config, for which there is no engine
// in this build yet — distinguished from a config error so a calling
// script can tell "your TOML is broken" from "your TOML is fine, this
// feature doesn't exist yet" apart.

#include <algorithm>
#include <array>
#include <chrono>
#include <cstdint>
#include <cstdlib>
#include <exception>
#include <filesystem>
#include <fstream>
#include <iomanip>
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
#include "engines/geometry_analysis/pore_analysis.hpp"
#include "engines/monte_carlo/framework_energy_grid.hpp"
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
        << "  aleator pore analyze <config.toml> [flags]      Analyze pore geometry (LCD, PLD, "
           "ASA, AV)\n"
        << "  aleator validate <config.toml>                  Validate a config, then exit\n"
        << "  aleator bench [--json]                          Run a quick built-in benchmark\n"
        << "\n"
        << "Only engines that are actually implemented get a subcommand: GCMC and\n"
        << "pore-geometry analysis. Molecular dynamics is not implemented yet, so there is\n"
        << "no `md run` -- see CLAUDE.md section 0 for status. `aleator validate` still\n"
        << "recognizes an [md] config file (the schema exists) and will say so plainly,\n"
        << "distinct from a malformed config.\n"
        << "\n"
        << "Flags (gcmc run, pore analyze):\n"
        << "  --dry-run   Validate the config and structure file, print the fully resolved\n"
        << "              configuration (every default included), and exit without\n"
        << "              running anything or writing any output.\n"
        << "  --json      Emit the final result as one line of JSON on stdout instead of\n"
        << "              human-readable text. Progress/diagnostic messages still go to\n"
        << "              stderr, so this is safe to pipe.\n"
        << "\n"
        << "A real (non-dry-run) `gcmc run` or `pore analyze` writes the fully resolved\n"
        << "configuration, including the RNG seed, to\n"
        << "<run.output_directory>/resolved_config.json -- enough to reproduce a published\n"
        << "result from its artifacts alone.\n"
        << "\n"
        << "Exit codes: 0 success, 1 configuration/usage error, 2 `validate` was pointed\n"
        << "at a valid config for an engine this build doesn't implement yet.\n";
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
    std::uint32_t adsorbateSpeciesIndex;
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

    return {std::move(structure),   std::move(forceField), std::move(adsorbate),
            adsorbateSpeciesIndex, fugacityPascal,         idealGas};
}

/// CLAUDE.md invariant #11 ("configuration keys are honored or rejected"):
/// `gcmc.energy_grid_spacing_angstrom`, if set, must actually build and use
/// a FrameworkEnergyGrid, not just be parsed and ignored. Deliberately not
/// called by prepareGcmc()/on --dry-run: construction is real,
/// O(gridPoints * frameworkAtomCount) work (see FrameworkEnergyGrid's doc
/// comment) that can itself take far longer than a validation pass should
/// (measured on the CH4/IRMOF-1 system: ~165 s at the spacing accurate
/// enough to matter -- CLAUDE.md section 5's "Resolved this session"), so
/// --dry-run intentionally does not pay it.
std::shared_ptr<const aleator::engines::FrameworkEnergyGrid> buildFrameworkEnergyGridIfConfigured(
    const aleator::io::GcmcConfig& gcmc, const PreparedGcmc& prepared) {
    if (!gcmc.energyGridSpacingAngstrom.has_value()) {
        return nullptr;
    }
    spdlog::info("building guest-host energy grid at {} Ang spacing (gcmc.energy_grid_spacing_angstrom "
                 "is set) -- this is one-time O(gridPoints * frameworkAtomCount) work...",
                 *gcmc.energyGridSpacingAngstrom);
    const auto buildStart = std::chrono::steady_clock::now();
    auto grid = std::make_shared<const aleator::engines::FrameworkEnergyGrid>(
        *prepared.forceField, prepared.structure.particles, prepared.structure.lattice,
        std::vector<std::uint32_t>{prepared.adsorbateSpeciesIndex}, *gcmc.energyGridSpacingAngstrom);
    const auto buildEnd = std::chrono::steady_clock::now();
    spdlog::info("energy grid built in {:.1f} s",
                 std::chrono::duration<double>(buildEnd - buildStart).count());
    return grid;
}

/// Single source of truth for the resolved-config JSON schema, shared by
/// `gcmc run --dry-run --json`'s stdout output and the real-run
/// `resolved_config.json` artifact (writeResolvedConfig) -- one serializer,
/// not two that could silently drift apart.
std::string resolvedConfigJson(const aleator::io::GcmcRunConfig& cfg) {
    const auto& run = cfg.run;
    const auto& g = cfg.gcmc;
    std::ostringstream oss;
    // Round-trip-safe precision (17 significant decimal digits is always
    // enough to recover the exact original IEEE754 double -- the default
    // stream precision of 6 is not, e.g. it silently truncates
    // mass_amu=16.04246 to 16.0425): this file is meant to let "a
    // published result be reproduced from its artifacts alone" (CLAUDE.md
    // CLI milestone), and truncated physics parameters (temperature,
    // pressure, cutoff, epsilon/sigma) would make that reproduction only
    // approximate.
    oss << std::setprecision(17);
    oss << "{\"run\":{\"name\":\"" << jsonEscape(run.name) << "\",\"output_directory\":\""
        << jsonEscape(run.outputDirectory.string()) << "\",\"rng_seed\":" << run.rngSeed
        << ",\"thread_count\":" << run.threadCount << "},"
        << "\"gcmc\":{\"framework_cif\":\"" << jsonEscape(g.frameworkCif.string())
        << "\",\"temperature_kelvin\":" << g.temperatureKelvin
        << ",\"pressure_bar\":" << g.pressureBar << ",\"cutoff_angstrom\":" << g.cutoffAngstrom
        << ",\"equilibration_steps\":" << g.equilibrationSteps
        << ",\"production_steps\":" << g.productionSteps << ",\"energy_grid_spacing_angstrom\":"
        << (g.energyGridSpacingAngstrom.has_value() ? std::to_string(*g.energyGridSpacingAngstrom)
                                                      : "null")
        << ",\"adsorbate\":{\"name\":\"" << jsonEscape(g.adsorbate.name)
        << "\",\"epsilon_kelvin\":" << g.adsorbate.epsilonKelvin
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
        oss << "{\"element\":\"" << jsonEscape(s.element) << "\",\"epsilon_kelvin\":" << s.epsilonKelvin
            << ",\"sigma_angstrom\":" << s.sigmaAngstrom << "}";
    }
    oss << "]}}";
    return oss.str();
}

void printGcmcConfig(const aleator::io::GcmcRunConfig& cfg, bool json) {
    const auto& run = cfg.run;
    const auto& g = cfg.gcmc;
    if (json) {
        std::cout << resolvedConfigJson(cfg) << "\n";
        return;
    }

    std::cout << std::setprecision(17);
    std::cout << "run.name                       = " << run.name << "\n"
              << "run.output_directory           = " << run.outputDirectory.string() << "\n"
              << "run.rng_seed                   = " << run.rngSeed << "\n"
              << "run.thread_count                = " << run.threadCount << "\n"
              << "gcmc.framework_cif             = " << g.frameworkCif.string() << "\n"
              << "gcmc.temperature_kelvin        = " << g.temperatureKelvin << "\n"
              << "gcmc.pressure_bar              = " << g.pressureBar << "\n"
              << "gcmc.cutoff_angstrom           = " << g.cutoffAngstrom << "\n"
              << "gcmc.equilibration_steps       = " << g.equilibrationSteps << "\n"
              << "gcmc.production_steps          = " << g.productionSteps << "\n"
              << "gcmc.energy_grid_spacing_angstrom = "
              << (g.energyGridSpacingAngstrom.has_value() ? std::to_string(*g.energyGridSpacingAngstrom)
                                                             : "not set (direct O(N) guest-host scan)")
              << "\n"
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

// CODATA 2018 / SI 2019 exact -- same constants and conversion
// tests/validation/test_gcmc_ch4_irmof1_isotherm.cc uses, so the CLI's
// printed mmol/g figure is directly comparable to that test's (and to
// published isotherms, which are conventionally reported in mmol/g, not
// molecules/Ang^3).
constexpr double kAmuToKg = 1.66053906660e-27;
constexpr double kAvogadro = 6.02214076e23;

/// Writes the fully-resolved configuration (as JSON, including
/// `run.rng_seed`) to `<outputDirectory>/resolved_config.json`, creating
/// the directory if needed. This is the artifact CLAUDE.md's milestone
/// asks for: "a published result can be reproduced from its artifacts
/// alone." Only called for a real run -- see printUsage()/runGcmcCommand's
/// doc comments for why --dry-run doesn't write anything.
void writeResolvedConfig(const aleator::io::GcmcRunConfig& cfg) {
    std::filesystem::create_directories(cfg.run.outputDirectory);
    const auto path = cfg.run.outputDirectory / "resolved_config.json";
    std::ofstream out(path);
    if (!out) {
        throw std::runtime_error("could not open " + path.string() + " for writing (run.output_directory "
                                  "= " +
                                  cfg.run.outputDirectory.string() + ")");
    }
    out << resolvedConfigJson(cfg);
    spdlog::info("wrote resolved configuration (incl. run.rng_seed={}) to {}", cfg.run.rngSeed,
                 path.string());
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

    const auto frameworkEnergyGrid = buildFrameworkEnergyGridIfConfigured(cfg.gcmc, prepared);

    // Computed before particles_ is moved into the engine below.
    double frameworkMassAmu = 0.0;
    for (double m : prepared.structure.particles.mass) {
        frameworkMassAmu += m;
    }
    const double frameworkMassKg = frameworkMassAmu * kAmuToKg;
    const double volume = std::abs(prepared.structure.lattice.volume());

    writeResolvedConfig(cfg);

    auto rng = std::make_unique<aleator::core::Philox4x32Rng>();
    rng->seed(cfg.run.rngSeed, /*streamIndex=*/0);
    aleator::engines::MonteCarloEngine engine(
        std::move(prepared.structure.particles), prepared.structure.lattice, prepared.forceField,
        std::move(rng), cfg.gcmc.temperatureKelvin, prepared.adsorbate, prepared.fugacityPascal,
        /*maxTranslationDisplacementAngstrom=*/1.5, /*maxRotationAngleRadians=*/0.5,
        /*electrostatics=*/nullptr, frameworkEnergyGrid);

    spdlog::info("equilibrating ({} steps)...", cfg.gcmc.equilibrationSteps);
    engine.run(cfg.gcmc.equilibrationSteps);

    spdlog::info("production ({} steps)...", cfg.gcmc.productionSteps);
    constexpr int kProgressUpdates = 10;
    const std::size_t stepsPerUpdate =
        std::max<std::size_t>(1, cfg.gcmc.productionSteps / kProgressUpdates);
    // Samples moleculeCount() after EVERY production step (not once per
    // progress-update batch): this is the same estimator
    // tests/validation/test_gcmc_ch4_irmof1_isotherm.cc and
    // tests/known_deviation/test_gcmc_ch4_irmof1_known_deviation.cc use
    // (their block-averaging is a standard-error computation over this
    // same per-step sample stream, not a coarser mean) -- deliberately, so
    // this command's printed <N> is directly comparable to, and with a
    // matching seed/config reproduces exactly, those tests' numbers
    // (CLAUDE.md milestone: "its CLI output must match the validation
    // test's numbers").
    double sumN = 0.0;
    std::size_t stepsRun = 0;
    while (stepsRun < cfg.gcmc.productionSteps) {
        const std::size_t thisBatch = std::min(stepsPerUpdate, cfg.gcmc.productionSteps - stepsRun);
        for (std::size_t i = 0; i < thisBatch; ++i) {
            engine.run(1);
            sumN += static_cast<double>(engine.moleculeCount());
        }
        stepsRun += thisBatch;
        spdlog::info("  [{:3.0f}%] step {}/{}  <N> so far = {:.3f}  molecules = {}",
                     100.0 * static_cast<double>(stepsRun) /
                         static_cast<double>(cfg.gcmc.productionSteps),
                     stepsRun, cfg.gcmc.productionSteps,
                     sumN / static_cast<double>(stepsRun), engine.moleculeCount());
    }
    const double meanN = sumN / static_cast<double>(cfg.gcmc.productionSteps);
    const double loadingMoleculesPerAngstrom3 = meanN / volume;
    // moles adsorbed per gram of framework, in mmol/g -- matching
    // tests/known_deviation/test_gcmc_ch4_irmof1_known_deviation.cc's
    // conversion exactly (kg->g and mol->mmol cancel, so it's just
    // (meanN / kAvogadro) / frameworkMassKg, not derived from the
    // molecules/Ang^3 density above).
    const double loadingMmolPerG =
        frameworkMassKg > 0.0 ? (meanN / kAvogadro) / frameworkMassKg : 0.0;

    if (json) {
        std::ostringstream oss;
        oss << "{\"mean_molecule_count\":" << meanN << ",\"loading_molecules_per_angstrom3\":"
            << loadingMoleculesPerAngstrom3 << ",\"loading_mmol_per_g\":" << loadingMmolPerG
            << ",\"cell_volume_angstrom3\":" << volume << "}";
        std::cout << oss.str() << "\n";
    } else {
        std::cout << "Final result:\n"
                  << "  <N> (mean molecule count over production) = " << meanN << "\n"
                  << "  loading                                   = " << loadingMoleculesPerAngstrom3
                  << " molecules/Ang^3\n"
                  << "                                             = " << loadingMmolPerG
                  << " mmol/g\n";
    }
    return kExitSuccess;
}

// ---------------------------------------------------------------- pore ---

// engines/geometry_analysis is now implemented (CLAUDE.md section 0), so
// `pore analyze` is a real subcommand -- `md run` still is not
// (engines/dynamics remains NotImplemented) and stays absent from --help.

/// Per-atom radii for `structure`, using this codebase's own Zeo++-sourced
/// default radius table (pore_analysis.hpp's
/// zeoPlusPlusDefaultRadiusAngstrom) -- deliberately not a config key: the
/// milestone that added this command reused the built-in table rather than
/// inventing a bespoke `[[pore.radii]]` schema for something the table
/// already covers for every element these example/test frameworks use.
/// Throws std::invalid_argument (via zeoPlusPlusDefaultRadiusAngstrom), by
/// element name, if the CIF contains an element the table doesn't have --
/// caught by main()'s generic handler like any other config-time error.
std::vector<double> poreRadiiFor(const aleator::io::StructureData& structure) {
    std::vector<double> radii;
    radii.reserve(structure.speciesSymbols.size());
    for (const auto& symbol : structure.speciesSymbols) {
        radii.push_back(aleator::engines::zeoPlusPlusDefaultRadiusAngstrom(symbol));
    }
    return radii;
}

std::string resolvedConfigJson(const aleator::io::PoreRunConfig& cfg) {
    const auto& run = cfg.run;
    const auto& p = cfg.pore;
    std::ostringstream oss;
    oss << std::setprecision(17);
    oss << "{\"run\":{\"name\":\"" << jsonEscape(run.name) << "\",\"output_directory\":\""
        << jsonEscape(run.outputDirectory.string()) << "\",\"rng_seed\":" << run.rngSeed
        << ",\"thread_count\":" << run.threadCount << "},"
        << "\"pore\":{\"framework_cif\":\"" << jsonEscape(p.frameworkCif.string())
        << "\",\"probe_radius_angstrom\":" << p.probeRadiusAngstrom << "}}";
    return oss.str();
}

void printPoreConfig(const aleator::io::PoreRunConfig& cfg, bool json) {
    if (json) {
        std::cout << resolvedConfigJson(cfg) << "\n";
        return;
    }
    const auto& run = cfg.run;
    const auto& p = cfg.pore;
    std::cout << std::setprecision(17);
    std::cout << "run.name                     = " << run.name << "\n"
              << "run.output_directory         = " << run.outputDirectory.string() << "\n"
              << "run.rng_seed                 = " << run.rngSeed << "\n"
              << "run.thread_count             = " << run.threadCount << "\n"
              << "pore.framework_cif           = " << p.frameworkCif.string() << "\n"
              << "pore.probe_radius_angstrom   = " << p.probeRadiusAngstrom << "\n";
}

/// Mirrors writeResolvedConfig(GcmcRunConfig) -- see that function's doc
/// comment for why this artifact exists.
void writeResolvedConfig(const aleator::io::PoreRunConfig& cfg) {
    std::filesystem::create_directories(cfg.run.outputDirectory);
    const auto path = cfg.run.outputDirectory / "resolved_config.json";
    std::ofstream out(path);
    if (!out) {
        throw std::runtime_error("could not open " + path.string() + " for writing (run.output_directory "
                                  "= " +
                                  cfg.run.outputDirectory.string() + ")");
    }
    out << resolvedConfigJson(cfg);
    spdlog::info("wrote resolved configuration (incl. run.rng_seed={}) to {}", cfg.run.rngSeed,
                 path.string());
}

int runPoreCommand(const std::filesystem::path& configPath, bool dryRun, bool json) {
    const aleator::io::PoreRunConfig cfg = aleator::io::loadPoreConfig(configPath);
    const auto structure = aleator::io::readCif(cfg.pore.frameworkCif);
    const std::vector<double> radii = poreRadiiFor(structure); // validates element coverage too

    if (dryRun) {
        printPoreConfig(cfg, json);
        return kExitSuccess;
    }

    spdlog::info("loaded {} framework atoms ({} species) from {}", structure.particles.size(),
                 structure.speciesSymbols.size(), cfg.pore.frameworkCif.string());

    aleator::engines::PoreAnalysisOptions options;
    options.probeRadiusAngstrom = cfg.pore.probeRadiusAngstrom;
    options.sampleSeed = cfg.run.rngSeed;

    spdlog::info("analyzing pore geometry (probe radius = {} Ang)...", options.probeRadiusAngstrom);
    const auto result =
        aleator::engines::analyzePoreGeometry(structure.particles, structure.lattice, radii, options);

    writeResolvedConfig(cfg);

    if (json) {
        std::ostringstream oss;
        oss << std::setprecision(17);
        oss << "{\"largest_cavity_diameter_angstrom\":" << result.largestCavityDiameterAngstrom
            << ",\"pore_limiting_diameter_angstrom\":" << result.poreLimitingDiameterAngstrom
            << ",\"accessible_surface_area_angstrom2\":" << result.accessibleSurfaceAreaAngstromSq
            << ",\"inaccessible_surface_area_angstrom2\":" << result.inaccessibleSurfaceAreaAngstromSq
            << ",\"accessible_volume_angstrom3\":" << result.accessibleVolumeAngstromCubed
            << ",\"inaccessible_volume_angstrom3\":" << result.inaccessibleVolumeAngstromCubed
            << ",\"total_volume_angstrom3\":" << result.totalVolumeAngstromCubed << "}";
        std::cout << oss.str() << "\n";
    } else {
        std::cout << "Final result:\n"
                  << "  largest cavity diameter (LCD)    = " << result.largestCavityDiameterAngstrom
                  << " Ang\n"
                  << "  pore limiting diameter (PLD)     = " << result.poreLimitingDiameterAngstrom
                  << " Ang\n"
                  << "  accessible surface area (ASA)    = " << result.accessibleSurfaceAreaAngstromSq
                  << " Ang^2\n"
                  << "  inaccessible surface area (NASA) = " << result.inaccessibleSurfaceAreaAngstromSq
                  << " Ang^2\n"
                  << "  accessible volume (AV)           = " << result.accessibleVolumeAngstromCubed
                  << " Ang^3\n"
                  << "  inaccessible volume (NAV)        = " << result.inaccessibleVolumeAngstromCubed
                  << " Ang^3\n"
                  << "  total unit cell volume           = " << result.totalVolumeAngstromCubed
                  << " Ang^3\n";
    }
    return kExitSuccess;
}

// ------------------------------------------------------------ validate ---

// `gcmc run` and `pore analyze` are real, executable subcommands
// (engines/geometry_analysis is now implemented); `md run` still is not
// (engines/dynamics remains NotImplemented). `validate` still recognizes
// an `[md]` config file -- its schema is real (io/config.hpp) and checking
// it is honest, useful groundwork, distinct from claiming the engine
// itself runs. Real schema validation happens first (so a genuinely
// malformed `[md]` config still gets a precise key/line error, same as
// `[gcmc]`/`[pore]`); only after that succeeds does this say plainly that
// there is no engine for it yet -- never silently "valid" in a way that
// could be mistaken for "runnable" (invariant #7: no stub presented as
// complete).
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
            poreRadiiFor(structure); // validates every element has a known radius too
            std::cout << configPath.string() << ": valid [pore] config (\"" << cfg.run.name << "\", "
                      << structure.particles.size() << " framework atoms)\n";
            return kExitSuccess;
        }
        case aleator::io::ConfigKind::Md: {
            const auto cfg = aleator::io::loadMdConfig(configPath);
            std::cout << configPath.string() << ": schema-valid [md] config (\"" << cfg.run.name
                      << "\"), but molecular dynamics is not implemented in this build yet -- there is "
                         "no `aleator md run`\n";
            return kExitNotImplemented;
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
                std::cerr << "error: usage: aleator pore analyze <config.toml> [--dry-run] [--json]\n";
                return kExitConfigError;
            }
            const auto parsed =
                parseFlags(std::vector<std::string_view>(rawArgs.begin() + 2, rawArgs.end()));
            if (parsed.positional.size() != 1) {
                std::cerr << "error: usage: aleator pore analyze <config.toml> [--dry-run] [--json]\n";
                return kExitConfigError;
            }
            return runPoreCommand(std::string(parsed.positional[0]), parsed.dryRun, parsed.json);
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
