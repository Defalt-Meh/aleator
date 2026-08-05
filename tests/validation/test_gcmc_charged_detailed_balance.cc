// Validation anchor for the charged-GCMC milestone: exercises the REAL
// MonteCarloEngine (not EwaldIncrementalState in isolation --
// tests/validation/test_ewald_incremental_state.cc already covers that)
// running actual translation/rotation/insertion/deletion moves on a
// charged system, via MonteCarloEngine::run(). Two things this test alone
// can catch that the isolated EwaldIncrementalState tests structurally
// cannot: a wiring bug in MonteCarloEngine's attempt*() methods (wrong
// site list, a forgotten commit call, propose/commit argument mismatch),
// and whether the real acceptance-criteria integration (dispersion +
// electrostatic deltaU combined) still round-trips correctly.
//
// CLAUDE.md milestone: "Run 10^5 GCMC steps with a charged species...
// compare the incrementally-maintained S(k) and total Ewald energy
// against a full from-scratch recomputation" -- done here through the
// actual production code path, complementing the isolated-state version.

#include <catch2/catch_test_macros.hpp>

#include <array>
#include <cmath>
#include <cstdio>
#include <memory>
#include <vector>

#include "core/geometry/lattice.hpp"
#include "core/math/counter_based_rng.hpp"
#include "core/math/particle_data.hpp"
#include "engines/monte_carlo/molecule_species.hpp"
#include "engines/monte_carlo/monte_carlo_engine.hpp"
#include "forcefield/electrostatics/ewald.hpp"
#include "forcefield/pairwise/lennard_jones.hpp"

using aleator::core::Lattice;
using aleator::core::ParticleData;
using aleator::core::Philox4x32Rng;
using aleator::engines::MoleculeSpecies;
using aleator::engines::MonteCarloEngine;
using aleator::forcefield::Ewald;
using aleator::forcefield::LennardJones;
using aleator::forcefield::LennardJonesParameters;
using aleator::forcefield::LennardJonesTruncation;

namespace {

Lattice cubicLattice(double length) {
    std::array<std::array<double, 3>, 3> m{};
    m[0] = {length, 0.0, 0.0};
    m[1] = {0.0, length, 0.0};
    m[2] = {0.0, 0.0, length};
    return Lattice(m);
}

/// Small charged framework, LJ species 0: alternating +1/-1 point charges,
/// net neutral -- same non-physical-but-sufficient construction as
/// test_ewald_incremental_state.cc's, reused here specifically to give the
/// REAL engine's dispersion force field something to act on too (a bare
/// point charge with zero LJ epsilon would let molecules collapse onto
/// the framework, which is fine for this test's purpose -- it only cares
/// about electrostatic cache correctness, not producing a physically
/// realistic isotherm).
ParticleData buildChargedFramework(double boxLength) {
    ParticleData framework;
    const std::array<double, 2> coords{boxLength * 0.25, boxLength * 0.75};
    int sign = 1;
    for (double x : coords) {
        for (double y : coords) {
            for (double z : coords) {
                framework.push_back(x, y, z, 0.0, 0.0, 0.0, /*mass=*/20.0,
                                     /*charge=*/static_cast<double>(sign), /*species=*/0);
                sign = -sign;
            }
        }
    }
    return framework;
}

} // namespace

TEST_CASE("Real charged GCMC (MonteCarloEngine, not EwaldIncrementalState in isolation) shows "
          "no cache drift over 10^5 steps",
          "[validation][montecarlo][electrostatics][slow]") {
    const double boxLength = 20.0;
    const Lattice lattice = cubicLattice(boxLength);
    ParticleData framework = buildChargedFramework(boxLength);

    // LJ species 0 = framework, 1 = CO2's C, 2 = CO2's O -- dispersion is
    // secondary to this test's purpose (electrostatic cache correctness),
    // so parameters are modest, not literature values.
    const std::vector<LennardJonesParameters> ljSpecies{
        {50.0, 3.0},  // 0: framework
        {30.0, 2.8},  // 1: C
        {80.0, 3.0},  // 2: O
    };
    const double ljCutoff = 8.0; // < L_perp/2 = 10.0 for this cubic box
    auto forceField =
        std::make_shared<const LennardJones>(ljSpecies, ljCutoff, LennardJonesTruncation::Shifted);

    const double realSpaceCutoff = 9.0; // < L_perp/2 = 10.0
    auto electrostatics = std::make_shared<const Ewald>(/*alpha=*/0.25, realSpaceCutoff,
                                                          /*kMax=*/6);

    MoleculeSpecies co2;
    co2.sites.push_back({0.0, 0.0, -1.16, /*ljSpecies=*/2, /*charge=*/-0.35, /*mass=*/15.9994});
    co2.sites.push_back({0.0, 0.0, 0.0, /*ljSpecies=*/1, /*charge=*/0.70, /*mass=*/12.0107});
    co2.sites.push_back({0.0, 0.0, 1.16, /*ljSpecies=*/2, /*charge=*/-0.35, /*mass=*/15.9994});

    auto rng = std::make_unique<Philox4x32Rng>();
    rng->seed(/*seed=*/2026, /*streamIndex=*/0);

    // A generous fugacity, chosen only to keep a nontrivial, bounded
    // population of molecules present throughout the run (this test cares
    // about cache correctness at realistic occupancy, not about producing
    // a physically meaningful loading).
    MonteCarloEngine engine(framework, lattice, forceField, std::move(rng),
                             /*temperatureKelvin=*/298.0, co2, /*fugacityPascal=*/5.0e6);

    constexpr std::size_t kTotalSteps = 100'000;
    constexpr std::size_t kChunk = 1'000;
    double maxDrift = 0.0;
    double maxRelDrift = 0.0;

    for (std::size_t done = 0; done < kTotalSteps; done += kChunk) {
        engine.run(kChunk);
        const double drift = engine.electrostaticEnergyDriftForTesting();
        const double energy = engine.electrostaticEnergyForTesting();
        const double relDrift = drift / std::max(1.0, std::abs(energy));
        maxDrift = std::max(maxDrift, drift);
        maxRelDrift = std::max(maxRelDrift, relDrift);
        INFO("after " << (done + kChunk) << " steps: drift=" << drift << " energy=" << energy
                       << " molecules=" << engine.moleculeCount());
        REQUIRE(relDrift <= 1e-10);
    }

    std::printf(
        "Real-engine charged GCMC drift gate: %zu steps (checked every %zu), max absolute "
        "drift = %.3e, max relative drift = %.3e, final molecule count = %zu\n",
        kTotalSteps, kChunk, maxDrift, maxRelDrift, engine.moleculeCount());
}
