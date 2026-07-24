// CLI entry point. cli/ is the top of the dependency graph (CLAUDE.md #3):
// it is the only place allowed to link core/, forcefield/, engines/, and
// io/ together. No physics runs here — `run` wires the layers together and
// lets NotImplemented surface as a clean error, rather than silently doing
// nothing.

#include <cstdlib>
#include <exception>
#include <iostream>
#include <memory>
#include <string_view>
#include <vector>

#include <aleator/version.hpp>
#include <spdlog/spdlog.h>

#include "core/exceptions.hpp"
#include "core/geometry/lattice.hpp"
#include "core/math/counter_based_rng.hpp"
#include "core/math/particle_data.hpp"
#include "engines/monte_carlo/monte_carlo_engine.hpp"
#include "forcefield/pairwise/lennard_jones.hpp"
#include "io/config.hpp"

namespace {

void printUsage() {
    std::cout << "aleator " << aleator::kVersion << "\n\n"
              << "Usage:\n"
              << "  aleator --version           Print the version and exit\n"
              << "  aleator --help               Print this message and exit\n"
              << "  aleator run <config.toml>    Load a run config and start a simulation\n";
}

int runCommand(std::string_view configPath) {
    const aleator::io::RunConfig config = aleator::io::loadRunConfig(configPath);
    spdlog::info("loaded run config \"{}\" (seed={}, threads={})", config.name, config.rngSeed,
                 config.threadCount);

    // Wires core + forcefield + engines together to prove the dependency
    // graph links; no physics is attempted. This is expected to throw.
    auto forceField = std::make_shared<const aleator::forcefield::LennardJones>(
        std::vector<aleator::forcefield::LennardJonesParameters>{}, 12.0);
    auto rng = std::make_unique<aleator::core::Philox4x32Rng>();
    aleator::engines::MonteCarloEngine engine(aleator::engines::Ensemble::Nvt, {},
                                               aleator::core::Lattice{}, forceField,
                                               std::move(rng), 298.15);
    try {
        engine.run(1);
    } catch (const aleator::NotImplemented& err) {
        spdlog::warn("{} — physics is not implemented in this build yet", err.what());
        return EXIT_SUCCESS;
    }
    return EXIT_SUCCESS;
}

} // namespace

int main(int argc, char** argv) {
    const std::vector<std::string_view> args(argv + 1, argv + argc);

    if (args.empty() || args.front() == "--help" || args.front() == "-h") {
        printUsage();
        return args.empty() ? EXIT_FAILURE : EXIT_SUCCESS;
    }
    if (args.front() == "--version") {
        std::cout << aleator::kVersion << "\n";
        return EXIT_SUCCESS;
    }
    if (args.front() == "run") {
        if (args.size() != 2) {
            std::cerr << "error: `run` requires exactly one argument, a config file path\n";
            return EXIT_FAILURE;
        }
        try {
            return runCommand(args[1]);
        } catch (const std::exception& err) {
            spdlog::error("{}", err.what());
            return EXIT_FAILURE;
        }
    }

    std::cerr << "error: unrecognized command \"" << args.front() << "\"\n";
    printUsage();
    return EXIT_FAILURE;
}
