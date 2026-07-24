// Integration test across core/ + forcefield/ + engines/: proves the
// dependency graph actually links and the objects can be constructed
// together end to end. No physics exists yet, so the only correct outcome
// is a clean NotImplemented — not a crash, not a silent success.

#include <catch2/catch_test_macros.hpp>

#include <memory>
#include <vector>

#include "core/exceptions.hpp"
#include "core/geometry/lattice.hpp"
#include "core/math/counter_based_rng.hpp"
#include "core/math/particle_data.hpp"
#include "engines/monte_carlo/monte_carlo_engine.hpp"
#include "forcefield/pairwise/lennard_jones.hpp"

TEST_CASE("MonteCarloEngine wires core + forcefield + engines and reports NotImplemented",
          "[integration][engines]") {
    auto forceField = std::make_shared<const aleator::forcefield::LennardJones>(
        std::vector<aleator::forcefield::LennardJonesParameters>{}, 12.0);
    auto rng = std::make_unique<aleator::core::Philox4x32Rng>();

    aleator::engines::MonteCarloEngine engine(aleator::engines::Ensemble::Nvt,
                                               aleator::core::ParticleData{},
                                               aleator::core::Lattice{}, forceField,
                                               std::move(rng), 298.15);

    REQUIRE(engine.ensemble() == aleator::engines::Ensemble::Nvt);
    REQUIRE_THROWS_AS(engine.run(1), aleator::NotImplemented);
}
