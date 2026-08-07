#pragma once

#include <array>
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
/// line/column, if `name` is missing, the file cannot be parsed as TOML, or
/// `thread_count` is set to anything other than 1 -- no threading is
/// implemented yet (CLAUDE.md invariant #11), so any other value would be
/// silently ignored rather than honored, and is rejected instead.
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

    /// CLAUDE.md section 5 performance milestone: grid spacing (Å) for the
    /// O(1)-interpolation guest-host energy field (engines/monte_carlo/
    /// framework_energy_grid.hpp), replacing the direct O(frameworkCount)
    /// scan. Unset (the default) means "don't build one -- use the direct
    /// scan", not "use some implicit spacing": this is a real accuracy/
    /// speed tradeoff (see FrameworkEnergyGrid's doc comment and
    /// tests/validation/test_framework_energy_grid.cc's measured
    /// spacing-vs-accuracy-vs-build-time sweep, and README.md), so it is
    /// opt-in rather than silently changed under an existing config
    /// (CLAUDE.md invariant #11: a key the user didn't set must not change
    /// what their run computes).
    std::optional<double> energyGridSpacingAngstrom;

    /// Directory FrameworkEnergyGrid cache files live in (see
    /// engines/monte_carlo/framework_energy_grid.hpp's
    /// loadOrBuildFrameworkEnergyGrid) -- only meaningful if
    /// energyGridSpacingAngstrom is set. Resolved relative to the config
    /// file's own directory (same convention as framework_cif/
    /// output_directory), default "energy_grid_cache" -- deliberately NOT
    /// inside run.output_directory by default, since a grid's whole point
    /// is amortization across pressure points/repeat runs/subsequent jobs
    /// (CLAUDE.md section 5 milestone), which usually means separate
    /// per-run output directories but a config file (or sibling configs)
    /// sharing the same framework.
    std::filesystem::path energyGridCacheDirectory = "energy_grid_cache";

    /// Optional override of the automatically-computed minimum supercell
    /// replication (core::minimumSupercellReplication, applied in
    /// src/cli/main.cc's prepareGcmc() once the framework CIF's real
    /// lattice is known -- this loader only checks each component is a
    /// positive integer; the CIF-dependent "is this at least the computed
    /// minimum" check happens there, once the lattice is available).
    /// Unset means "use the automatically-computed minimum" -- the
    /// override can only ever raise the replication above that computed
    /// minimum, never below it (a below-minimum override would violate
    /// minimum image, CLAUDE.md invariant #10, and is rejected).
    std::optional<std::array<int, 3>> supercellOverride;
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

/// Which of `[gcmc]` / `[pore]` a config file's top level contains — used by
/// `aleator validate` to dispatch to the right loader (and its full
/// validation) without the caller having to already know what kind of config
/// it's pointing at.
enum class ConfigKind { Gcmc, Pore, Unknown };

[[nodiscard]] ConfigKind detectConfigKind(const std::filesystem::path& file);

} // namespace aleator::io
