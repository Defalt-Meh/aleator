#include "io/config.hpp"

#include <sstream>
#include <string>
#include <utility>

#include <toml++/toml.hpp>

#include "io/config_error.hpp"

namespace aleator::io {

namespace {

/// (line, column) of `node`, or of `fallback` if `node` is null (the "key
/// doesn't exist at all" case, where toml++ has nothing more precise to
/// offer than the enclosing table's own position).
std::pair<std::uint32_t, std::uint32_t> locationOf(const toml::node* node, const toml::table& fallback) {
    const auto& src = node != nullptr ? node->source() : fallback.source();
    return {src.begin.line, src.begin.column};
}

std::string typeName(const toml::node& node) {
    std::ostringstream oss;
    oss << node.type();
    return oss.str();
}

const toml::table& requireTable(const toml::table& parent, std::string_view key,
                                 const std::string& filePath, const std::string& keyPath) {
    const auto view = parent[key];
    const auto* tbl = view.as_table();
    if (tbl == nullptr) {
        const auto [line, column] = locationOf(view.node(), parent);
        const std::string problem = view.node() == nullptr
                                         ? "missing required table"
                                         : "expected a table, found " + typeName(*view.node());
        throw ConfigError(filePath, keyPath, problem, line, column);
    }
    return *tbl;
}

template <typename T>
T requireValue(const toml::table& parent, std::string_view key, const std::string& filePath,
               const std::string& keyPath, std::string_view typeDescription) {
    const auto view = parent[key];
    const auto value = view.template value<T>();
    if (!value.has_value()) {
        const auto [line, column] = locationOf(view.node(), parent);
        const std::string problem =
            view.node() == nullptr
                ? "missing required key (expected " + std::string(typeDescription) + ")"
                : "expected " + std::string(typeDescription) + ", found " + typeName(*view.node());
        throw ConfigError(filePath, keyPath, problem, line, column);
    }
    return *value;
}

template <typename T>
T optionalValue(const toml::table& parent, std::string_view key, T defaultValue,
                 const std::string& filePath, const std::string& keyPath,
                 std::string_view typeDescription) {
    const auto view = parent[key];
    if (view.node() == nullptr) {
        return defaultValue;
    }
    const auto value = view.template value<T>();
    if (!value.has_value()) {
        const auto [line, column] = locationOf(view.node(), parent);
        throw ConfigError(filePath, keyPath,
                           "expected " + std::string(typeDescription) + ", found " +
                               typeName(*view.node()),
                           line, column);
    }
    return *value;
}

void requirePositive(double value, const std::string& filePath, const std::string& keyPath,
                      const toml::table& parent, std::string_view key) {
    if (!(value > 0.0)) {
        const auto [line, column] = locationOf(parent[key].node(), parent);
        throw ConfigError(filePath, keyPath, "must be positive (got " + std::to_string(value) + ")",
                           line, column);
    }
}

std::filesystem::path resolveRelativeTo(const std::filesystem::path& configFile,
                                         const std::string& relativeOrAbsolute) {
    const std::filesystem::path candidate(relativeOrAbsolute);
    if (candidate.is_absolute()) {
        return candidate.lexically_normal();
    }
    return (configFile.parent_path() / candidate).lexically_normal();
}

toml::table parseFile(const std::filesystem::path& file) {
    try {
        return toml::parse_file(file.string());
    } catch (const toml::parse_error& err) {
        const auto& pos = err.source().begin;
        throw ConfigError(file.string(), "<parse>", std::string(err.description()), pos.line,
                           pos.column);
    }
}

RunConfig parseRunTable(const toml::table& root, const std::filesystem::path& file) {
    const std::string filePath = file.string();
    const toml::table& run = requireTable(root, "run", filePath, "run");

    RunConfig config;
    config.name = requireValue<std::string>(run, "name", filePath, "run.name", "a string");
    // Resolved relative to the config file's own directory, same as
    // framework_cif/structure_file below -- not the process's current
    // working directory, so that "a published result can be reproduced
    // from its artifacts alone" (CLAUDE.md CLI milestone) doesn't also
    // silently depend on which directory `aleator` happened to be invoked
    // from.
    config.outputDirectory = resolveRelativeTo(
        file, optionalValue<std::string>(run, "output_directory", std::string("."), filePath,
                                          "run.output_directory", "a string"));
    config.rngSeed = static_cast<std::uint64_t>(
        optionalValue<std::int64_t>(run, "rng_seed", 0, filePath, "run.rng_seed", "an integer"));
    config.threadCount = static_cast<unsigned>(optionalValue<std::int64_t>(
        run, "thread_count", 1, filePath, "run.thread_count", "an integer"));
    // CLAUDE.md invariant #11: "A key the config schema accepts and the
    // code ignores is a lie to the user." There is no threading anywhere
    // in src/ (see CLAUDE.md section 0) -- every run is single-threaded
    // regardless of this value, so any value other than the one that
    // matches reality is rejected here rather than silently accepted and
    // ignored.
    if (config.threadCount != 1) {
        const auto [line, column] = locationOf(run["thread_count"].node(), run);
        throw ConfigError(filePath, "run.thread_count",
                           "not supported in this build: no threading is implemented yet "
                           "(every run is single-threaded regardless of this value); only "
                           "thread_count = 1 is accepted",
                           line, column);
    }
    return config;
}

} // namespace

RunConfig loadRunConfig(const std::filesystem::path& file) {
    const toml::table root = parseFile(file);
    return parseRunTable(root, file);
}

GcmcRunConfig loadGcmcConfig(const std::filesystem::path& file) {
    const toml::table root = parseFile(file);
    const std::string filePath = file.string();

    GcmcRunConfig result;
    result.run = parseRunTable(root, file);

    const toml::table& gcmc = requireTable(root, "gcmc", filePath, "gcmc");
    GcmcConfig& cfg = result.gcmc;

    const std::string cifPath =
        requireValue<std::string>(gcmc, "framework_cif", filePath, "gcmc.framework_cif", "a string");
    cfg.frameworkCif = resolveRelativeTo(file, cifPath);

    cfg.temperatureKelvin = requireValue<double>(gcmc, "temperature_kelvin", filePath,
                                                  "gcmc.temperature_kelvin", "a number");
    requirePositive(cfg.temperatureKelvin, filePath, "gcmc.temperature_kelvin", gcmc,
                     "temperature_kelvin");

    cfg.pressureBar =
        requireValue<double>(gcmc, "pressure_bar", filePath, "gcmc.pressure_bar", "a number");
    requirePositive(cfg.pressureBar, filePath, "gcmc.pressure_bar", gcmc, "pressure_bar");

    cfg.cutoffAngstrom = optionalValue<double>(gcmc, "cutoff_angstrom", cfg.cutoffAngstrom, filePath,
                                                "gcmc.cutoff_angstrom", "a number");
    requirePositive(cfg.cutoffAngstrom, filePath, "gcmc.cutoff_angstrom", gcmc, "cutoff_angstrom");

    const auto equilibrationSteps = optionalValue<std::int64_t>(
        gcmc, "equilibration_steps", static_cast<std::int64_t>(cfg.equilibrationSteps), filePath,
        "gcmc.equilibration_steps", "an integer");
    const auto productionSteps = optionalValue<std::int64_t>(
        gcmc, "production_steps", static_cast<std::int64_t>(cfg.productionSteps), filePath,
        "gcmc.production_steps", "an integer");
    if (equilibrationSteps < 0) {
        const auto [line, column] = locationOf(gcmc["equilibration_steps"].node(), gcmc);
        throw ConfigError(filePath, "gcmc.equilibration_steps", "must not be negative", line, column);
    }
    if (productionSteps <= 0) {
        const auto [line, column] = locationOf(gcmc["production_steps"].node(), gcmc);
        throw ConfigError(filePath, "gcmc.production_steps", "must be positive", line, column);
    }
    cfg.equilibrationSteps = static_cast<std::size_t>(equilibrationSteps);
    cfg.productionSteps = static_cast<std::size_t>(productionSteps);

    if (gcmc["energy_grid_spacing_angstrom"].node() != nullptr) {
        cfg.energyGridSpacingAngstrom = requireValue<double>(
            gcmc, "energy_grid_spacing_angstrom", filePath, "gcmc.energy_grid_spacing_angstrom",
            "a number");
        requirePositive(*cfg.energyGridSpacingAngstrom, filePath,
                         "gcmc.energy_grid_spacing_angstrom", gcmc, "energy_grid_spacing_angstrom");
    }
    cfg.energyGridCacheDirectory = resolveRelativeTo(
        file, optionalValue<std::string>(gcmc, "energy_grid_cache_directory",
                                          cfg.energyGridCacheDirectory.string(), filePath,
                                          "gcmc.energy_grid_cache_directory", "a string"));

    // Supercell replication override (CLAUDE.md milestone: "automatic
    // supercell replication... a config key to override the computed
    // replication upward, never downward"). Only syntax/positivity is
    // checked here -- whether this override is actually >= the computed
    // minimum requires the framework CIF's real lattice, which this
    // parse-only loader doesn't read (framework_lj's per-element coverage
    // check has the same structural reason for living in
    // src/cli/main.cc's prepareGcmc() instead of here).
    if (gcmc["supercell"].node() != nullptr) {
        const auto* supercellArray = gcmc["supercell"].as_array();
        if (supercellArray == nullptr || supercellArray->size() != 3) {
            const auto [line, column] = locationOf(gcmc["supercell"].node(), gcmc);
            throw ConfigError(filePath, "gcmc.supercell",
                               "expected an array of exactly 3 positive integers [nx, ny, nz]", line,
                               column);
        }
        std::array<int, 3> supercell{};
        for (std::size_t i = 0; i < 3; ++i) {
            const auto value = (*supercellArray)[i].value<std::int64_t>();
            if (!value.has_value() || *value < 1) {
                const auto [line, column] = locationOf(&(*supercellArray)[i], gcmc);
                throw ConfigError(filePath, "gcmc.supercell",
                                   "each entry must be a positive integer (>= 1)", line, column);
            }
            supercell[i] = static_cast<int>(*value);
        }
        cfg.supercellOverride = supercell;
    }

    const auto* frameworkLjArray = gcmc["framework_lj"].as_array();
    if (frameworkLjArray == nullptr || frameworkLjArray->empty()) {
        const auto [line, column] = locationOf(gcmc["framework_lj"].node(), gcmc);
        throw ConfigError(filePath, "gcmc.framework_lj",
                           "missing or empty; at least one [[gcmc.framework_lj]] entry is required",
                           line, column);
    }
    for (std::size_t i = 0; i < frameworkLjArray->size(); ++i) {
        const auto* entry = (*frameworkLjArray)[i].as_table();
        const std::string entryPath = "gcmc.framework_lj[" + std::to_string(i) + "]";
        if (entry == nullptr) {
            const auto& src = (*frameworkLjArray)[i].source();
            throw ConfigError(filePath, entryPath, "expected a table", src.begin.line,
                               src.begin.column);
        }
        GcmcFrameworkSpecies species;
        species.element = requireValue<std::string>(*entry, "element", filePath,
                                                      entryPath + ".element", "a string");
        species.epsilonKelvin = requireValue<double>(*entry, "epsilon_kelvin", filePath,
                                                       entryPath + ".epsilon_kelvin", "a number");
        requirePositive(species.epsilonKelvin, filePath, entryPath + ".epsilon_kelvin", *entry,
                         "epsilon_kelvin");
        species.sigmaAngstrom = requireValue<double>(*entry, "sigma_angstrom", filePath,
                                                       entryPath + ".sigma_angstrom", "a number");
        requirePositive(species.sigmaAngstrom, filePath, entryPath + ".sigma_angstrom", *entry,
                         "sigma_angstrom");
        cfg.frameworkLennardJones.push_back(std::move(species));
    }

    const toml::table& adsorbate = requireTable(gcmc, "adsorbate", filePath, "gcmc.adsorbate");
    cfg.adsorbate.name =
        requireValue<std::string>(adsorbate, "name", filePath, "gcmc.adsorbate.name", "a string");
    cfg.adsorbate.epsilonKelvin = requireValue<double>(
        adsorbate, "epsilon_kelvin", filePath, "gcmc.adsorbate.epsilon_kelvin", "a number");
    requirePositive(cfg.adsorbate.epsilonKelvin, filePath, "gcmc.adsorbate.epsilon_kelvin", adsorbate,
                     "epsilon_kelvin");
    cfg.adsorbate.sigmaAngstrom = requireValue<double>(
        adsorbate, "sigma_angstrom", filePath, "gcmc.adsorbate.sigma_angstrom", "a number");
    requirePositive(cfg.adsorbate.sigmaAngstrom, filePath, "gcmc.adsorbate.sigma_angstrom", adsorbate,
                     "sigma_angstrom");
    cfg.adsorbate.massAmu =
        requireValue<double>(adsorbate, "mass_amu", filePath, "gcmc.adsorbate.mass_amu", "a number");
    requirePositive(cfg.adsorbate.massAmu, filePath, "gcmc.adsorbate.mass_amu", adsorbate, "mass_amu");
    cfg.adsorbate.chargeE = optionalValue<double>(adsorbate, "charge_e", 0.0, filePath,
                                                   "gcmc.adsorbate.charge_e", "a number");

    // Critical properties are all-or-nothing: partially specifying them is
    // almost certainly a typo (a missing key silently falling back to
    // "use ideal gas" while the other two are honored would be a confusing
    // half-state), so it's rejected as a config error rather than guessed.
    const bool hasTc = adsorbate["critical_temperature_kelvin"].node() != nullptr;
    const bool hasPc = adsorbate["critical_pressure_pascal"].node() != nullptr;
    const bool hasOmega = adsorbate["acentric_factor"].node() != nullptr;
    if (hasTc || hasPc || hasOmega) {
        if (!(hasTc && hasPc && hasOmega)) {
            const auto [line, column] = locationOf(&adsorbate, gcmc);
            throw ConfigError(filePath, "gcmc.adsorbate",
                               "critical_temperature_kelvin, critical_pressure_pascal, and "
                               "acentric_factor must all be given together, or all omitted "
                               "(omitted -> ideal-gas fugacity approximation)",
                               line, column);
        }
        GcmcCriticalProperties props;
        props.criticalTemperatureKelvin =
            requireValue<double>(adsorbate, "critical_temperature_kelvin", filePath,
                                  "gcmc.adsorbate.critical_temperature_kelvin", "a number");
        requirePositive(props.criticalTemperatureKelvin, filePath,
                         "gcmc.adsorbate.critical_temperature_kelvin", adsorbate,
                         "critical_temperature_kelvin");
        props.criticalPressurePascal =
            requireValue<double>(adsorbate, "critical_pressure_pascal", filePath,
                                  "gcmc.adsorbate.critical_pressure_pascal", "a number");
        requirePositive(props.criticalPressurePascal, filePath,
                         "gcmc.adsorbate.critical_pressure_pascal", adsorbate,
                         "critical_pressure_pascal");
        props.acentricFactor = requireValue<double>(
            adsorbate, "acentric_factor", filePath, "gcmc.adsorbate.acentric_factor", "a number");
        cfg.adsorbate.criticalProperties = props;
    }

    return result;
}

PoreRunConfig loadPoreConfig(const std::filesystem::path& file) {
    const toml::table root = parseFile(file);
    const std::string filePath = file.string();

    PoreRunConfig result;
    result.run = parseRunTable(root, file);

    const toml::table& pore = requireTable(root, "pore", filePath, "pore");
    const std::string cifPath =
        requireValue<std::string>(pore, "framework_cif", filePath, "pore.framework_cif", "a string");
    result.pore.frameworkCif = resolveRelativeTo(file, cifPath);
    result.pore.probeRadiusAngstrom =
        optionalValue<double>(pore, "probe_radius_angstrom", result.pore.probeRadiusAngstrom,
                               filePath, "pore.probe_radius_angstrom", "a number");
    requirePositive(result.pore.probeRadiusAngstrom, filePath, "pore.probe_radius_angstrom", pore,
                     "probe_radius_angstrom");

    return result;
}

ConfigKind detectConfigKind(const std::filesystem::path& file) {
    const toml::table root = parseFile(file);
    const bool hasGcmc = root["gcmc"].as_table() != nullptr;
    const bool hasPore = root["pore"].as_table() != nullptr;
    if (hasGcmc && !hasPore) {
        return ConfigKind::Gcmc;
    }
    if (hasPore && !hasGcmc) {
        return ConfigKind::Pore;
    }
    return ConfigKind::Unknown;
}

} // namespace aleator::io
