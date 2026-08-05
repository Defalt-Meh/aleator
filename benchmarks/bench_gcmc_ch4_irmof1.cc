// CLAUDE.md section 5 performance milestone: tracked baseline for the GCMC
// hot loop (MonteCarloEngine::run), at the real CH4/IRMOF-1 validation
// system's scale (424 framework atoms, 25.832 Ang cubic cell, 12.0 Ang
// cutoff -- tests/validation/data/irmof1/). Uses a synthetic framework at
// that same atom count/density/cutoff rather than reading the real CIF
// (keeping this target self-contained, matching bench_neighbor_list.cc's
// convention), since what's measured here (per-step cost of
// MonteCarloEngine::run(), including this milestone's guest-guest CellList
// wiring and Lattice::minimumImageDisplacement's orthorhombic fast path) is
// governed by the same geometry either way.
//
// Before this milestone: profiled directly (not benchmarked, since this
// benchmark target and the guest-guest CellList/orthorhombic fast path
// didn't exist yet) at ~91% of total time inside
// Lattice::minimumImageDisplacement's 125-candidate triclinic search,
// called from a direct O(frameworkCount) scan for every trial move -- see
// README.md's performance section. The real end-to-end number this
// benchmark's workload corresponds to: the actual 4-point CH4/IRMOF-1
// known-deviation isotherm test
// (tests/known_deviation/test_gcmc_ch4_irmof1_known_deviation.cc, 480,000
// total MC steps across 4 pressure points) went from 32.0s to 3.97s
// (Release, this machine) after this milestone's guest-guest CellList
// wiring plus the orthorhombic minimum-image fast path -- the latter, not
// the former, is what actually closed most of the gap: the CellList wiring
// alone made no measurable difference on this system (guest occupancy here
// tops out around 23 molecules, too few for an O(N) vs. O(1)-ish
// distinction to show up against a 424-atom framework scan).

#include <array>
#include <memory>
#include <random>
#include <vector>

#include <benchmark/benchmark.h>

#include "core/geometry/lattice.hpp"
#include "core/math/counter_based_rng.hpp"
#include "core/math/particle_data.hpp"
#include "engines/monte_carlo/molecule_species.hpp"
#include "engines/monte_carlo/monte_carlo_engine.hpp"
#include "forcefield/pairwise/lennard_jones.hpp"

namespace {

using aleator::core::CounterBasedRng;
using aleator::core::Lattice;
using aleator::core::ParticleData;
using aleator::core::Philox4x32Rng;
using aleator::engines::MoleculeSpecies;
using aleator::engines::MonteCarloEngine;
using aleator::forcefield::LennardJones;
using aleator::forcefield::LennardJonesParameters;
using aleator::forcefield::LennardJonesTruncation;

constexpr std::size_t kFrameworkAtomCount = 424; // matches real IRMOF-1.cif
constexpr double kCellLengthAngstrom = 25.832;    // matches real IRMOF-1.cif
constexpr double kCutoffAngstrom = 12.0;          // matches the validation tests

ParticleData randomFramework(const Lattice& lattice, unsigned seed) {
    std::mt19937 rng(seed);
    std::uniform_real_distribution<double> unit(0.0, 1.0);
    ParticleData framework;
    framework.resize(kFrameworkAtomCount);
    for (std::size_t i = 0; i < kFrameworkAtomCount; ++i) {
        const auto cart =
            lattice.fractionalToCartesian({unit(rng), unit(rng), unit(rng)});
        framework.x[i] = cart[0];
        framework.y[i] = cart[1];
        framework.z[i] = cart[2];
        framework.mass[i] = 65.38; // representative (Zn-scale) framework atom mass
        framework.charge[i] = 0.0;
        framework.species[i] = 0;
    }
    return framework;
}

void BM_GcmcRun(benchmark::State& state) {
    std::array<std::array<double, 3>, 3> m{};
    m[0] = {kCellLengthAngstrom, 0.0, 0.0};
    m[1] = {0.0, kCellLengthAngstrom, 0.0};
    m[2] = {0.0, 0.0, kCellLengthAngstrom};
    const Lattice lattice(m);

    const std::vector<LennardJonesParameters> ljParameters{
        {62.4, 2.46},  // framework species, UFF-Zn-scale
        {158.5, 3.72}, // guest species, CH4_sp3-scale
    };
    auto forceField = std::make_shared<const LennardJones>(ljParameters, kCutoffAngstrom,
                                                             LennardJonesTruncation::Shifted);

    MoleculeSpecies guest;
    guest.sites.push_back({0.0, 0.0, 0.0, /*ljSpecies=*/1, /*charge=*/0.0, /*mass=*/16.04});

    const auto numSteps = static_cast<std::size_t>(state.range(0));

    for (auto _ : state) {
        state.PauseTiming();
        auto framework = randomFramework(lattice, 42);
        auto rng = std::make_unique<Philox4x32Rng>();
        rng->seed(/*seed=*/7, /*streamIndex=*/0);
        MonteCarloEngine engine(std::move(framework), lattice, forceField, std::move(rng),
                                 /*temperatureKelvin=*/298.0, guest,
                                 /*fugacityPascal=*/1.0e5);
        state.ResumeTiming();

        engine.run(numSteps);
        benchmark::DoNotOptimize(engine.moleculeCount());
    }
    state.SetItemsProcessed(static_cast<std::int64_t>(state.iterations()) *
                             static_cast<std::int64_t>(numSteps));
}

} // namespace

// 120,000 steps matches one pressure point's equilibration + production in
// the real known-deviation isotherm test; 480,000 matches all four points
// combined.
BENCHMARK(BM_GcmcRun)->Arg(120'000)->Arg(480'000)->Unit(benchmark::kMillisecond)->Iterations(1);
