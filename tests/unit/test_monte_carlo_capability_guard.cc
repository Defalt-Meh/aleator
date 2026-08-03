// CLAUDE.md invariant #10: "Interfaces that a consumer calls are
// implemented by every concrete class the consumer can be handed... a
// base-class default that throws NotImplemented is a runtime hole
// disguised as a compile-time contract... the engine must reject loudly
// [at construction]." MonteCarloEngine's trial moves all call
// ForceField::computeParticleEnergy(), which Ewald does not implement
// (see forcefield/electrostatics/ewald.hpp) -- this test proves the
// engine's constructor catches that mismatch immediately, by name, rather
// than the caller finding out only when a move throws NotImplemented
// mid-run.

#include <catch2/catch_test_macros.hpp>

#include <array>
#include <memory>
#include <string>

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

namespace {

Lattice cubicLattice(double length) {
    std::array<std::array<double, 3>, 3> m{};
    m[0] = {length, 0.0, 0.0};
    m[1] = {0.0, length, 0.0};
    m[2] = {0.0, 0.0, length};
    return Lattice(m);
}

} // namespace

TEST_CASE("ForceField capability defaults and overrides are self-consistent", "[unit][forcefield]") {
    const std::vector<LennardJonesParameters> species{{100.0, 3.0}};
    const LennardJones lj(species, 6.0);
    REQUIRE(lj.name() == "LennardJones");
    REQUIRE(lj.supportsSingleParticleEnergy());

    const Ewald ewald(/*alpha=*/0.3, /*realSpaceCutoff=*/8.0, /*kMax=*/6);
    REQUIRE(ewald.name() == "Ewald");
    REQUIRE_FALSE(ewald.supportsSingleParticleEnergy());
}

TEST_CASE("MonteCarloEngine's GCMC constructor rejects a force field lacking "
          "single-particle-energy support, by name",
          "[unit][montecarlo]") {
    const double boxLength = 20.0;
    const Lattice lattice = cubicLattice(boxLength);
    ParticleData framework;
    framework.push_back(1.0, 1.0, 1.0, 1.0, 0.0, 0.0, /*mass=*/22.99, /*charge=*/1.0,
                         /*species=*/0);
    framework.push_back(2.0, 2.0, 2.0, -1.0, 0.0, 0.0, /*mass=*/35.45, /*charge=*/-1.0,
                         /*species=*/0);

    auto forceField = std::make_shared<const Ewald>(/*alpha=*/0.3, /*realSpaceCutoff=*/8.0,
                                                      /*kMax=*/6);

    MoleculeSpecies adsorbate;
    adsorbate.sites.push_back({0.0, 0.0, 0.0, /*ljSpecies=*/0, /*charge=*/0.0, /*mass=*/16.0});

    auto rng = std::make_unique<Philox4x32Rng>();
    rng->seed(/*seed=*/1, /*streamIndex=*/0);

    try {
        MonteCarloEngine engine(framework, lattice, forceField, std::move(rng),
                                 /*temperatureKelvin=*/298.0, adsorbate,
                                 /*fugacityPascal=*/1.0e5);
        FAIL("expected MonteCarloEngine construction to throw for an Ewald force field");
    } catch (const std::invalid_argument& err) {
        const std::string what = err.what();
        REQUIRE(what.find("Ewald") != std::string::npos);
        REQUIRE(what.find("supportsSingleParticleEnergy") != std::string::npos);
    }
}
