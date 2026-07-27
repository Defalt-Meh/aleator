// forcefield/pairwise milestone: mixing rules and truncation-scheme
// mechanics that don't need an external reference to check — hand-computed
// values and the defining continuity properties of each scheme. The NIST
// SRSW energy/virial validation lives in
// tests/validation/test_lennard_jones_nist.cc.

#include <catch2/catch_test_macros.hpp>

#include <cmath>
#include <stdexcept>
#include <vector>

#include "forcefield/pairwise/lennard_jones.hpp"

using aleator::forcefield::LennardJones;
using aleator::forcefield::LennardJonesParameters;
using aleator::forcefield::LennardJonesTruncation;
using aleator::forcefield::MixingRule;

namespace {
constexpr double kEps = 1e-12;
}

TEST_CASE("mixedParameters returns the exact per-species value for i == j", "[unit][forcefield]") {
    std::vector<LennardJonesParameters> species{{1.5, 2.0}, {3.0, 4.0}};
    LennardJones lj(species, 10.0);
    const auto mixed = lj.mixedParameters(1, 1);
    REQUIRE(mixed.epsilon == 3.0);
    REQUIRE(mixed.sigma == 4.0);
}

TEST_CASE("Lorentz-Berthelot mixing: arithmetic mean sigma, geometric mean epsilon",
          "[unit][forcefield]") {
    std::vector<LennardJonesParameters> species{{1.0, 2.0}, {4.0, 6.0}};
    LennardJones lj(species, 10.0, LennardJonesTruncation::Truncated, MixingRule::LorentzBerthelot);
    const auto mixed = lj.mixedParameters(0, 1);
    REQUIRE(std::abs(mixed.epsilon - std::sqrt(1.0 * 4.0)) < kEps); // = 2.0
    REQUIRE(std::abs(mixed.sigma - 0.5 * (2.0 + 6.0)) < kEps);      // = 4.0
}

TEST_CASE("Geometric mixing: geometric mean for both sigma and epsilon", "[unit][forcefield]") {
    std::vector<LennardJonesParameters> species{{1.0, 2.0}, {4.0, 8.0}};
    LennardJones lj(species, 10.0, LennardJonesTruncation::Truncated, MixingRule::Geometric);
    const auto mixed = lj.mixedParameters(0, 1);
    REQUIRE(std::abs(mixed.epsilon - std::sqrt(1.0 * 4.0)) < kEps); // = 2.0
    REQUIRE(std::abs(mixed.sigma - std::sqrt(2.0 * 8.0)) < kEps);   // = 4.0
}

TEST_CASE("mixedParameters throws on an out-of-range species id", "[unit][forcefield]") {
    std::vector<LennardJonesParameters> species{{1.0, 1.0}};
    LennardJones lj(species, 10.0);
    REQUIRE_THROWS_AS(lj.mixedParameters(0, 5), std::out_of_range);
    REQUIRE_THROWS_AS(lj.mixedParameters(5, 0), std::out_of_range);
}

TEST_CASE("LennardJones constructor rejects a non-positive cutoff", "[unit][forcefield]") {
    std::vector<LennardJonesParameters> species{{1.0, 1.0}};
    REQUIRE_THROWS_AS(LennardJones(species, 0.0), std::invalid_argument);
    REQUIRE_THROWS_AS(LennardJones(species, -1.0), std::invalid_argument);
}

TEST_CASE("Truncated: energy and force are exactly zero beyond the cutoff", "[unit][forcefield]") {
    std::vector<LennardJonesParameters> species{{1.0, 1.0}};
    LennardJones lj(species, 3.0, LennardJonesTruncation::Truncated);
    const auto& p = species[0];
    REQUIRE(lj.pairEnergy(3.0001, p) == 0.0);
    REQUIRE(lj.pairForceOverR(3.0001, p) == 0.0);
    REQUIRE(lj.pairEnergy(2.9999, p) != 0.0);
}

TEST_CASE("Truncated and Shifted have identical forces (a constant energy offset doesn't change "
          "the derivative)",
          "[unit][forcefield]") {
    std::vector<LennardJonesParameters> species{{1.0, 1.0}};
    LennardJones truncated(species, 3.0, LennardJonesTruncation::Truncated);
    LennardJones shifted(species, 3.0, LennardJonesTruncation::Shifted);
    const auto& p = species[0];
    for (double r : {0.9, 1.0, 1.5, 2.0, 2.9}) {
        REQUIRE(std::abs(truncated.pairForceOverR(r, p) - shifted.pairForceOverR(r, p)) < kEps);
    }
}

TEST_CASE("Shifted: energy is continuous (zero) at the cutoff", "[unit][forcefield]") {
    std::vector<LennardJonesParameters> species{{1.0, 1.0}};
    LennardJones lj(species, 3.0, LennardJonesTruncation::Shifted);
    const auto& p = species[0];
    // Just inside the cutoff, the shifted energy should approach 0.
    REQUIRE(std::abs(lj.pairEnergy(2.999999, p)) < 1e-5);
}

TEST_CASE("LinearForceShifted: both energy and force are continuous (zero) at the cutoff",
          "[unit][forcefield]") {
    std::vector<LennardJonesParameters> species{{1.0, 1.0}};
    LennardJones lj(species, 3.0, LennardJonesTruncation::LinearForceShifted);
    const auto& p = species[0];
    REQUIRE(std::abs(lj.pairEnergy(2.999999, p)) < 1e-5);
    REQUIRE(std::abs(lj.pairForceOverR(2.999999, p)) < 1e-4);
}
