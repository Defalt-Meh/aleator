// Validation anchor for the mixture-GCMC milestone's item 3:
//
//   "EWALD CORRECTNESS ACROSS SPECIES. Insertion of species A and deletion
//   of species B change the charge set differently. Run the existing drift
//   gate with a charged mixture."
//
// Same drift-gate methodology as
// tests/validation/test_gcmc_charged_detailed_balance.cc (which, despite
// its filename, is the real-engine drift-gate integration test, not the
// acceptance-ratio algebra check -- see that file's own header comment),
// extended to TWO charged species with genuinely different site counts,
// charge magnitudes, and masses: CO2 (3-site, charges +0.70/-0.35/-0.35)
// and a simple charged diatomic (2-site, charges +0.482/-0.482, not a real
// N2 model -- N2 itself has no permanent charges; this is a deliberately
// distinct synthetic charge distribution, not a physically-labelled
// species). Insertion/deletion/swap of either species now perturbs
// EwaldIncrementalState's committed cache with a DIFFERENT-shaped charge
// delta than the single-species test ever exercised (3 charges appearing/
// disappearing vs. 2, at different magnitudes) -- exactly the scenario
// item 3 asks to be checked, not merely asserted correct by the general
// computeDelta()/proposeMove() argument (that argument is about what
// EwaldIncrementalState's API can accept in principle; this test is about
// whether MonteCarloEngine actually wires it correctly for a mixture in
// practice).

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

/// Identical construction to test_gcmc_charged_detailed_balance.cc's
/// buildChargedFramework -- reused verbatim (not physically meaningful,
/// just a real net-neutral, non-trivial charge distribution for the real
/// dispersion force field and Ewald sum to act on).
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

TEST_CASE("Real charged MIXTURE GCMC (two species with different charge sets) shows no cache "
          "drift over 10^5 steps",
          "[validation][montecarlo][electrostatics][mixture][slow]") {
    const double boxLength = 20.0;
    const Lattice lattice = cubicLattice(boxLength);
    ParticleData framework = buildChargedFramework(boxLength);

    // LJ species 0 = framework, 1/2 = CO2's C/O, 3/4 = the diatomic's two
    // sites -- dispersion parameters are modest, not literature values
    // (this test's purpose is electrostatic cache correctness under a
    // mixture, not a physically meaningful isotherm).
    const std::vector<LennardJonesParameters> ljSpecies{
        {50.0, 3.0},  // 0: framework
        {30.0, 2.8},  // 1: CO2 C
        {80.0, 3.0},  // 2: CO2 O
        {60.0, 3.3},  // 3: diatomic site a
        {60.0, 3.3},  // 4: diatomic site b
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

    // A different site count (2, not 3), different charge magnitude
    // (0.482, not 0.70/0.35), and different mass -- deliberately NOT a
    // scaled copy of CO2, so a bug that only manifests for a
    // differently-shaped charge delta has a real chance to show up.
    MoleculeSpecies diatomic;
    diatomic.sites.push_back({0.0, 0.0, -0.55, /*ljSpecies=*/3, /*charge=*/0.482, /*mass=*/14.0});
    diatomic.sites.push_back({0.0, 0.0, 0.55, /*ljSpecies=*/4, /*charge=*/-0.482, /*mass=*/14.0});

    auto rng = std::make_unique<Philox4x32Rng>();
    rng->seed(/*seed=*/20260812, /*streamIndex=*/0);

    // Generous fugacities, chosen only to keep a nontrivial, bounded
    // population of BOTH species present throughout the run -- this test
    // cares about cache correctness at realistic occupancy, not about
    // producing a physically meaningful loading or selectivity.
    MonteCarloEngine engine(framework, lattice, forceField, std::move(rng),
                             /*temperatureKelvin=*/298.0, std::vector<MoleculeSpecies>{co2, diatomic},
                             std::vector<double>{5.0e6, 5.0e6}, /*maxTranslation=*/1.5,
                             /*maxRotation=*/0.5, electrostatics);

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
                       << " N_co2=" << engine.moleculeCountOfSpecies(0)
                       << " N_diatomic=" << engine.moleculeCountOfSpecies(1));
        REQUIRE(relDrift <= 1e-10);
    }

    std::printf("Real-engine charged MIXTURE GCMC drift gate: %zu steps (checked every %zu), max "
                "absolute drift = %.3e, max relative drift = %.3e, final N_co2=%zu, "
                "N_diatomic=%zu\n",
                kTotalSteps, kChunk, maxDrift, maxRelDrift, engine.moleculeCountOfSpecies(0),
                engine.moleculeCountOfSpecies(1));
}
