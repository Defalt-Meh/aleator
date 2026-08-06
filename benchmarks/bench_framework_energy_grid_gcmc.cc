// CLAUDE.md section 5's FrameworkEnergyGrid-serialization milestone: "BENCHMARK, and report
// all four: single cell, no grid (the current 3.97 s number) / single cell, grid (expect it to
// lose -- confirm) / 2x2x2 supercell, no grid / 2x2x2 supercell, grid, cold cache and warm
// cache reported separately."
//
// "single cell, no grid" is bench_gcmc_ch4_irmof1.cc's own BM_GcmcRun -- not duplicated here.
// This file adds the other three (four, counting cold/warm separately), all on the REAL
// IRMOF-1.cif system (not the synthetic framework bench_gcmc_ch4_irmof1.cc uses), because the
// grid-accuracy sweep this milestone also ran (tests/validation/test_framework_energy_grid.cc)
// needs the real force field/cutoff/spacing to mean anything, and reusing the same real system
// here keeps the grid-build-cost numbers comparable to that sweep's own measurements.
//
// Grid construction (or cache load) is deliberately excluded from each BENCHMARK's own timed
// region (state.PauseTiming/ResumeTiming, same convention bench_gcmc_ch4_irmof1.cc already
// uses for framework/engine construction) and reported separately via std::printf instead --
// the two costs answer different questions ("what does one grid build/load cost" vs. "what
// does the GCMC loop cost once amortized") and conflating them into one number would hide
// exactly the amortization tradeoff this milestone exists to characterize.

#include <array>
#include <chrono>
#include <cstdint>
#include <cstdio>
#include <filesystem>
#include <map>
#include <memory>
#include <stdexcept>
#include <string>
#include <vector>

#include <benchmark/benchmark.h>

#include "core/geometry/lattice.hpp"
#include "core/geometry/supercell.hpp"
#include "core/math/counter_based_rng.hpp"
#include "core/math/particle_data.hpp"
#include "engines/monte_carlo/framework_energy_grid.hpp"
#include "engines/monte_carlo/molecule_species.hpp"
#include "engines/monte_carlo/monte_carlo_engine.hpp"
#include "forcefield/pairwise/lennard_jones.hpp"
#include "io/structure_io.hpp"

namespace {

using aleator::core::Lattice;
using aleator::core::minimumSupercellReplication;
using aleator::core::ParticleData;
using aleator::core::Philox4x32Rng;
using aleator::core::replicateSupercell;
using aleator::engines::FrameworkEnergyGrid;
using aleator::engines::FrameworkEnergyGridCacheRequest;
using aleator::engines::loadOrBuildFrameworkEnergyGrid;
using aleator::engines::MoleculeSpecies;
using aleator::engines::MonteCarloEngine;
using aleator::forcefield::LennardJones;
using aleator::forcefield::LennardJonesParameters;
using aleator::forcefield::LennardJonesTruncation;

// Same accuracy-validated spacing test_framework_energy_grid.cc's sweep settled on for both
// the single unit cell AND the 2x2x2 supercell (this milestone's own re-run of that sweep on
// the supercell found the needed spacing UNCHANGED at fixed accuracy, not coarser -- see
// CLAUDE.md section 0 and README.md's performance section for that finding in full).
constexpr double kGridSpacingAngstrom = 0.2;
constexpr double kCutoffAngstrom = 12.0;

std::filesystem::path irmof1DataFile(const std::string& name) {
    return std::filesystem::path(ALEATOR_IRMOF1_DATA_DIR) / name;
}

LennardJonesParameters uffParameters(const std::string& element) {
    static const std::map<std::string, LennardJonesParameters> table{
        {"Zn", {62.3992, 2.46155}},
        {"O", {48.1581, 3.03315}},
        {"C", {47.8562, 3.47299}},
        {"H", {7.64893, 2.84642}},
    };
    const auto it = table.find(element);
    if (it == table.end()) {
        throw std::runtime_error("uffParameters: no GenericMOFs LJ parameters for element '" +
                                  element + "'");
    }
    return it->second;
}

struct System {
    ParticleData framework;
    Lattice lattice;
    std::shared_ptr<const LennardJones> forceField;
    std::uint32_t ch4SpeciesIndex = 0;
};

System buildSystem(int nx, int ny, int nz) {
    const auto unitCell = aleator::io::readCif(irmof1DataFile("IRMOF-1.cif"));

    System system;
    if (nx == 1 && ny == 1 && nz == 1) {
        system.framework = unitCell.particles;
        system.lattice = unitCell.lattice;
    } else {
        auto replicated = replicateSupercell(unitCell.particles, unitCell.lattice, nx, ny, nz);
        system.framework = std::move(replicated.particles);
        system.lattice = replicated.lattice;
    }

    std::vector<LennardJonesParameters> ljParameters;
    for (const auto& symbol : unitCell.speciesSymbols) {
        ljParameters.push_back(uffParameters(symbol));
    }
    system.ch4SpeciesIndex = static_cast<std::uint32_t>(ljParameters.size());
    ljParameters.push_back({158.5, 3.72}); // CH4_sp3, RASPA2 GenericMOFs

    system.forceField = std::make_shared<const LennardJones>(ljParameters, kCutoffAngstrom,
                                                               LennardJonesTruncation::Shifted);
    return system;
}

void runEngine(const System& system, std::shared_ptr<const FrameworkEnergyGrid> grid,
               std::size_t numSteps) {
    MoleculeSpecies guest;
    guest.sites.push_back(
        {0.0, 0.0, 0.0, system.ch4SpeciesIndex, /*charge=*/0.0, /*mass=*/16.04246});

    auto rng = std::make_unique<Philox4x32Rng>();
    rng->seed(/*seed=*/7, /*streamIndex=*/0);
    MonteCarloEngine engine(system.framework, system.lattice, system.forceField, std::move(rng),
                             /*temperatureKelvin=*/298.0, guest, /*fugacityPascal=*/1.0e5,
                             /*maxTranslation=*/1.5, /*maxRotation=*/0.5, /*electrostatics=*/nullptr,
                             std::move(grid));
    engine.run(numSteps);
    benchmark::DoNotOptimize(engine.moleculeCount());
}

// Scenario: single unit cell, grid. "Expect it to lose" per the milestone -- this benchmark
// exists to confirm that, with a fresh, current measurement, not the prior session's number.
void BM_GcmcRunSingleCellGrid(benchmark::State& state) {
    const System system = buildSystem(1, 1, 1);
    const auto numSteps = static_cast<std::size_t>(state.range(0));

    bool printedBuildTime = false;
    for (auto _ : state) {
        state.PauseTiming();
        const auto buildStart = std::chrono::steady_clock::now();
        auto grid = std::make_shared<const FrameworkEnergyGrid>(
            *system.forceField, system.framework, system.lattice,
            std::vector<std::uint32_t>{system.ch4SpeciesIndex}, kGridSpacingAngstrom);
        const auto buildEnd = std::chrono::steady_clock::now();
        if (!printedBuildTime) {
            std::printf("[single-cell grid] build time: %.3f s\n",
                        std::chrono::duration<double>(buildEnd - buildStart).count());
            printedBuildTime = true;
        }
        state.ResumeTiming();

        runEngine(system, grid, numSteps);
    }
    state.SetItemsProcessed(static_cast<std::int64_t>(state.iterations()) *
                             static_cast<std::int64_t>(numSteps));
}

// Scenario: 2x2x2 supercell, no grid -- the direct O(frameworkCount) scan at the much larger
// atom count a supercell forces cutoff << L into (see README.md's performance section for why
// this system exists as its own scenario rather than reusing bench_gcmc_ch4_irmof1.cc's
// synthetic single-cell system).
void BM_GcmcRunSupercellNoGrid(benchmark::State& state) {
    const System system = buildSystem(2, 2, 2);
    const auto numSteps = static_cast<std::size_t>(state.range(0));

    for (auto _ : state) {
        runEngine(system, /*grid=*/nullptr, numSteps);
    }
    state.SetItemsProcessed(static_cast<std::int64_t>(state.iterations()) *
                             static_cast<std::int64_t>(numSteps));
}

// Scenario: 2x2x2 supercell, grid, cold cache then warm cache -- the actual disk-cache code
// path (loadOrBuildFrameworkEnergyGrid), not a bare in-memory construction, so this measures
// exactly what a real `aleator gcmc run` pays on a cache miss vs. a cache hit.
void BM_GcmcRunSupercellGridColdThenWarm(benchmark::State& state) {
    const System system = buildSystem(2, 2, 2);
    const auto numSteps = static_cast<std::size_t>(state.range(0));

    const std::filesystem::path cacheDir =
        std::filesystem::temp_directory_path() / "aleator_bench_supercell_grid_cache";
    std::filesystem::remove_all(cacheDir); // guarantee a real cold start, not a stale leftover

    const FrameworkEnergyGridCacheRequest request{
        cacheDir, {system.ch4SpeciesIndex}, kGridSpacingAngstrom};

    bool printedTimings = false;
    for (auto _ : state) {
        state.PauseTiming();
        bool coldHit = false;
        const auto coldStart = std::chrono::steady_clock::now();
        auto coldGrid = std::make_shared<const FrameworkEnergyGrid>(
            loadOrBuildFrameworkEnergyGrid(*system.forceField, system.framework, system.lattice,
                                            request, &coldHit));
        const auto coldEnd = std::chrono::steady_clock::now();

        bool warmHit = false;
        const auto warmStart = std::chrono::steady_clock::now();
        auto warmGrid = std::make_shared<const FrameworkEnergyGrid>(
            loadOrBuildFrameworkEnergyGrid(*system.forceField, system.framework, system.lattice,
                                            request, &warmHit));
        const auto warmEnd = std::chrono::steady_clock::now();

        if (!printedTimings) {
            std::printf("[2x2x2 supercell grid] cold (%s): %.3f s, warm (%s): %.3f s\n",
                        coldHit ? "HIT -- unexpected" : "miss, built+cached",
                        std::chrono::duration<double>(coldEnd - coldStart).count(),
                        warmHit ? "hit, loaded from disk" : "MISS -- unexpected",
                        std::chrono::duration<double>(warmEnd - warmStart).count());
            printedTimings = true;
        }
        state.ResumeTiming();

        runEngine(system, warmGrid, numSteps);
    }
    state.SetItemsProcessed(static_cast<std::int64_t>(state.iterations()) *
                             static_cast<std::int64_t>(numSteps));

    std::filesystem::remove_all(cacheDir);
}

} // namespace

// Step counts deliberately smaller than bench_gcmc_ch4_irmof1.cc's 480,000 -- the supercell's
// per-step direct-scan cost is ~8x higher (proportional to framework atom count), and the
// 0.2 Ang grid build/load costs (printed separately, not part of these Arg() timings) already
// dominate wall-clock time for the grid scenarios. 20,000 steps is enough to get a stable
// per-step estimate without an impractically long benchmark run. Reported alongside
// SetItemsProcessed so ms/step is directly comparable across all four scenarios despite the
// differing step counts and differing framework sizes.
BENCHMARK(BM_GcmcRunSingleCellGrid)->Arg(20'000)->Unit(benchmark::kMillisecond)->Iterations(1);
BENCHMARK(BM_GcmcRunSupercellNoGrid)->Arg(20'000)->Unit(benchmark::kMillisecond)->Iterations(1);
BENCHMARK(BM_GcmcRunSupercellGridColdThenWarm)
    ->Arg(20'000)
    ->Unit(benchmark::kMillisecond)
    ->Iterations(1);
