#pragma once

#include <cstdint>
#include <filesystem>
#include <optional>
#include <string>

namespace aleator::io {

/// The subset of run configuration that is pure bookkeeping (run name,
/// output directory, RNG seed, thread count) rather than physics. Wired up
/// now, for real, to prove the TOML dependency and the general shape of
/// config loading; simulation-specific sections (ensemble parameters, force
/// field selection, ...) get added as those features land.
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
/// back to RunConfig's defaults. Throws std::runtime_error if `name` is
/// missing or the file cannot be parsed as TOML.
[[nodiscard]] RunConfig loadRunConfig(const std::filesystem::path& file);

} // namespace aleator::io
