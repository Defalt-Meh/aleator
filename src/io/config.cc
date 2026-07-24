#include "io/config.hpp"

#include <stdexcept>

#include <toml++/toml.hpp>

namespace aleator::io {

RunConfig loadRunConfig(const std::filesystem::path& file) {
    toml::table table;
    try {
        table = toml::parse_file(file.string());
    } catch (const toml::parse_error& err) {
        throw std::runtime_error("failed to parse " + file.string() + ": " +
                                  std::string(err.description()));
    }

    const auto* run = table["run"].as_table();
    if (run == nullptr) {
        throw std::runtime_error(file.string() + ": missing required [run] table");
    }

    RunConfig config;

    const auto name = (*run)["name"].value<std::string>();
    if (!name.has_value()) {
        throw std::runtime_error(file.string() + ": [run].name is required");
    }
    config.name = *name;

    if (const auto outputDirectory = (*run)["output_directory"].value<std::string>();
        outputDirectory.has_value()) {
        config.outputDirectory = *outputDirectory;
    }
    if (const auto rngSeed = (*run)["rng_seed"].value<std::int64_t>(); rngSeed.has_value()) {
        config.rngSeed = static_cast<std::uint64_t>(*rngSeed);
    }
    if (const auto threadCount = (*run)["thread_count"].value<std::int64_t>();
        threadCount.has_value()) {
        config.threadCount = static_cast<unsigned>(*threadCount);
    }

    return config;
}

} // namespace aleator::io
