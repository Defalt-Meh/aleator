// Integration test across io/ + the vendored TOML parser + the filesystem —
// not a physics test, but a real round trip (write a file, parse it back)
// rather than a mock.

#include <catch2/catch_test_macros.hpp>

#include <filesystem>
#include <fstream>
#include <stdexcept>

#include "io/config.hpp"

TEST_CASE("loadRunConfig round-trips a TOML file written to disk", "[integration][io]") {
    const auto path = std::filesystem::temp_directory_path() / "aleator_test_run_config.toml";
    {
        std::ofstream out(path);
        out << "[run]\n"
            << "name = \"smoke-test\"\n"
            << "output_directory = \"out\"\n"
            << "rng_seed = 12345\n"
            << "thread_count = 4\n";
    }

    const aleator::io::RunConfig config = aleator::io::loadRunConfig(path);

    REQUIRE(config.name == "smoke-test");
    REQUIRE(config.outputDirectory == "out");
    REQUIRE(config.rngSeed == 12345);
    REQUIRE(config.threadCount == 4);

    std::filesystem::remove(path);
}

TEST_CASE("loadRunConfig rejects a file with no [run] table", "[integration][io]") {
    const auto path = std::filesystem::temp_directory_path() / "aleator_test_bad_config.toml";
    {
        std::ofstream out(path);
        out << "[not_run]\nname = \"x\"\n";
    }

    REQUIRE_THROWS_AS(aleator::io::loadRunConfig(path), std::runtime_error);

    std::filesystem::remove(path);
}
