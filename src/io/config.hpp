#pragma once

#include <cstdint>
#include <filesystem>
#include <optional>
#include <string>
#include <vector>

namespace aleator::io {

/// The subset of run configuration that is pure bookkeeping (run name,
/// output directory, RNG seed, thread count) rather than physics — common
/// to every subcommand, always read from the top-level `[run]` table.
struct RunConfig {
    std::string name;
    std::filesystem::path outputDirectory = ".";
    std::uint64_t rngSeed = 0;
    unsigned threadCount = 1;
};

/// Reads a RunConfig from a TOML file with a top-level `[run]` table:
///
///     [run]
///     name = "example"
///     output_directory = "out"
///     rng_seed = 42
///     thread_count = 4
///
/// `output_directory`, `rng_seed`, and `thread_count` are optional and fall
/// back to RunConfig's defaults. Throws ConfigError (a std::runtime_error)
/// naming the offending key and, where toml++ can supply one, the exact
/// line/column, if `name` is missing or the file cannot be parsed as TOML.
[[nodiscard]] RunConfig loadRunConfig(const std::filesystem::path& file);

/// One framework element's Lennard-Jones parameters, from a
/// `[[gcmc.framework_lj]]` array-of-tables entry. Deliberately explicit in
/// every example config (CLAUDE.md's "no bespoke input language, no hidden
/// magic" spirit) rather than looked up from a built-in force-field
/// database this milestone doesn't add.
struct GcmcFrameworkSpecies {
    std::string element;
    double epsilonKelvin;
    double sigmaAngstrom;
};

/// Peng-Robinson critical properties for computing this adsorbate's
/// fugacity from the configured pressure (see engines/monte_carlo/
/// peng_robinson.hpp). Optional: if omitted, the CLI falls back to the
/// ideal-gas approximation (fugacity = pressure) rather than silently
/// assuming methane's constants for whatever species was actually
/// configured — logged clearly either way, never silent.
struct GcmcCriticalProperties {
    double criticalTemperatureKelvin;
    double criticalPressurePascal;
    double acentricFactor;
};

/// The single (currently: exactly one) adsorbate species GCMC inserts and
/// deletes, from `[gcmc.adsorbate]`. Single-site only at the CLI level
/// (matching the one production validation target, united-atom methane);
/// the engine itself supports multi-site rigid molecules, but wiring that
/// through the config schema is future scope.
struct GcmcAdsorbate {
    std::string name;
    double epsilonKelvin;
    double sigmaAngstrom;
    double massAmu;
    double chargeE = 0.0;
    std::optional<GcmcCriticalProperties> criticalProperties;
};

/// GCMC run parameters, from `[gcmc]` (+ its `adsorbate` subtable and
/// `framework_lj` array of tables).
struct GcmcConfig {
    std::filesystem::path frameworkCif;
    double temperatureKelvin;
    double pressureBar;
    double cutoffAngstrom = 12.0;
    std::size_t equilibrationSteps = 10'000;
    std::size_t productionSteps = 20'000;
    std::vector<GcmcFrameworkSpecies> frameworkLennardJones;
    GcmcAdsorbate adsorbate;
};

struct GcmcRunConfig {
    RunConfig run;
    GcmcConfig gcmc;
};

/// Reads and fully validates a GCMC run config (see `[[gcmc.framework_lj]]`
/// / `[gcmc.adsorbate]` above and `examples/gcmc_ch4_irmof1.toml` for the
/// documented schema by example). Every required key is checked and every
/// value range-checked (positive temperature, positive pressure, a cutoff
/// that isn't absurd, at least one framework species, step counts > 0)
/// before this function returns — the whole point being that a malformed
/// config fails here, at startup, with a precise error, rather than
/// mid-run. Throws ConfigError.
[[nodiscard]] GcmcRunConfig loadGcmcConfig(const std::filesystem::path& file);

/// Pore-geometry analysis parameters, from `[pore]`. The analysis engine
/// itself is not implemented yet (engines/geometry_analysis) — this loader
/// exists so `aleator pore analyze --dry-run` and config validation are
/// already real and useful ahead of that.
struct PoreConfig {
    std::filesystem::path frameworkCif;
    double probeRadiusAngstrom = 1.2;
};

struct PoreRunConfig {
    RunConfig run;
    PoreConfig pore;
};

[[nodiscard]] PoreRunConfig loadPoreConfig(const std::filesystem::path& file);

/// Molecular dynamics run parameters, from `[md]`. The integrator itself is
/// not implemented yet (engines/dynamics) — same rationale as PoreConfig.
struct MdConfig {
    std::filesystem::path structureFile;
    double timestepPicoseconds;
    std::size_t numSteps;
};

struct MdRunConfig {
    RunConfig run;
    MdConfig md;
};

[[nodiscard]] MdRunConfig loadMdConfig(const std::filesystem::path& file);

/// Which of `[gcmc]` / `[pore]` / `[md]` a config file's top level
/// contains — used by `aleator validate` to dispatch to the right loader
/// (and its full validation) without the caller having to already know
/// what kind of config it's pointing at.
enum class ConfigKind { Gcmc, Pore, Md, Unknown };

[[nodiscard]] ConfigKind detectConfigKind(const std::filesystem::path& file);

} // namespace aleator::io
