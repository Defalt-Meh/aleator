// Validation anchor (CLAUDE.md milestone: "Tail correction consistency
// between energy and pressure"). NIST's published reference table gives
// U_LRC but not a separate virial/pressure tail correction, so there's no
// external number to check tailVirialCorrection() against directly — the
// ground truth here is independent numerical quadrature of the same
// defining integrals (Allen & Tildesley eq. 2.98/2.99-style), computed
// completely separately from the closed-form formulas in lennard_jones.cc.
//
// Both integrals (energy: r^2 V(r); virial: r^3 * (-dV/dr)) are improper
// (rc to infinity). Substituting u = 1/r turns each into a plain
// polynomial in u over the finite, well-behaved interval [0, 1/rc] — no
// truncation-at-a-large-radius approximation needed, so Simpson's rule
// with enough subdivisions is essentially exact.

#include <catch2/catch_test_macros.hpp>

#include <array>
#include <cmath>
#include <numbers>
#include <utility>
#include <vector>

#include "core/geometry/lattice.hpp"
#include "core/math/particle_data.hpp"
#include "forcefield/pairwise/lennard_jones.hpp"

namespace {

using aleator::core::Lattice;
using aleator::core::ParticleData;
using aleator::forcefield::LennardJones;
using aleator::forcefield::LennardJonesParameters;

Lattice cubicLattice(double length) {
    std::array<std::array<double, 3>, 3> m{};
    m[0] = {length, 0.0, 0.0};
    m[1] = {0.0, length, 0.0};
    m[2] = {0.0, 0.0, length};
    return Lattice(m);
}

/// Composite Simpson's rule on [a, b] with `intervals` (must be even)
/// subdivisions.
template <typename F>
double simpson(F&& f, double a, double b, int intervals) {
    const double h = (b - a) / static_cast<double>(intervals);
    double sum = f(a) + f(b);
    for (int k = 1; k < intervals; ++k) {
        const double x = a + static_cast<double>(k) * h;
        sum += (k % 2 == 0 ? 2.0 : 4.0) * f(x);
    }
    return sum * h / 3.0;
}

/// Independent (not reusing any LennardJones code) numerical evaluation of
/// the per-pair tail integrals via the u = 1/r substitution:
///   energyIntegral = integral_rc^inf r^2 V(r) dr
///   virialIntegral = integral_rc^inf r^3 * (-dV/dr) dr
struct TailIntegrals {
    double energyIntegral;
    double virialIntegral;
};

TailIntegrals referenceTailIntegrals(double epsilon, double sigma, double rc) {
    const double sigma6 = std::pow(sigma, 6);
    const double sigma12 = sigma6 * sigma6;

    // V(1/u)/u^4 = 4*eps*sigma12*u^8 - 4*eps*sigma6*u^2
    auto energyIntegrand = [&](double u) {
        return 4.0 * epsilon * sigma12 * std::pow(u, 8) - 4.0 * epsilon * sigma6 * u * u;
    };
    // -dV/dr(1/u) / u^5 = -24*eps*sigma6*u^2 + 48*eps*sigma12*u^8
    auto virialIntegrand = [&](double u) {
        return -24.0 * epsilon * sigma6 * u * u + 48.0 * epsilon * sigma12 * std::pow(u, 8);
    };

    constexpr int kIntervals = 200000;
    const double uMax = 1.0 / rc;
    return {
        simpson(energyIntegrand, 0.0, uMax, kIntervals),
        simpson(virialIntegrand, 0.0, uMax, kIntervals),
    };
}

/// U_LRC_total / W_LRC_total for a single species from the raw integrals,
/// computed independently of LennardJones::tailEnergyCorrection /
/// tailVirialCorrection: (1/2) * N * 4*pi*rho * integral.
struct ReferenceTailCorrection {
    double energy;
    double virial;
};

ReferenceTailCorrection referenceTailCorrection(double epsilon, double sigma, double rc,
                                                 std::size_t n, double volume) {
    const auto integrals = referenceTailIntegrals(epsilon, sigma, rc);
    const double rho = static_cast<double>(n) / volume;
    const double prefactor = 0.5 * static_cast<double>(n) * 4.0 * std::numbers::pi * rho;
    return {prefactor * integrals.energyIntegral, prefactor * integrals.virialIntegral};
}

} // namespace

TEST_CASE("LennardJones tail energy/virial corrections match independent numerical quadrature",
          "[validation][forcefield]") {
    struct Case {
        double epsilon;
        double sigma;
        double rc;
        std::size_t n;
        double boxLength;
    };
    const std::vector<Case> cases{
        {1.0, 1.0, 3.0, 800, 10.0}, {1.0, 1.0, 4.0, 400, 10.0}, {1.0, 1.0, 2.5, 200, 8.0},
        {0.7, 1.3, 3.2, 500, 12.0},
    };

    for (const auto& c : cases) {
        INFO("epsilon=" << c.epsilon << " sigma=" << c.sigma << " rc=" << c.rc << " N=" << c.n);
        const std::vector<LennardJonesParameters> species{{c.epsilon, c.sigma}};
        const LennardJones lj(species, c.rc);

        ParticleData particles;
        particles.resize(c.n);
        for (auto& s : particles.species) {
            s = 0;
        }
        const Lattice lattice = cubicLattice(c.boxLength);
        const double volume = c.boxLength * c.boxLength * c.boxLength;

        const double productionEnergy = lj.tailEnergyCorrection(particles, lattice);
        const double productionVirial = lj.tailVirialCorrection(particles, lattice);
        const auto reference = referenceTailCorrection(c.epsilon, c.sigma, c.rc, c.n, volume);

        const double energyRelError =
            std::abs(productionEnergy - reference.energy) / std::abs(reference.energy);
        const double virialRelError =
            std::abs(productionVirial - reference.virial) / std::abs(reference.virial);
        REQUIRE(energyRelError < 1e-9);
        REQUIRE(virialRelError < 1e-9);
    }
}

TEST_CASE("LennardJones tail corrections generalize correctly to multiple species",
          "[validation][forcefield]") {
    // Two species with different (epsilon, sigma); reference computed by
    // summing the independent per-species-pair quadrature directly,
    // exercising the same N_a*N_b cross-term structure as the production
    // multi-species formula but via a completely separate code path.
    const std::vector<LennardJonesParameters> species{{1.0, 1.0}, {0.5, 1.5}};
    const double rc = 3.0;
    const LennardJones lj(species, rc);

    const std::size_t nA = 300;
    const std::size_t nB = 200;
    ParticleData particles;
    particles.resize(nA + nB);
    for (std::size_t i = 0; i < nA; ++i) {
        particles.species[i] = 0;
    }
    for (std::size_t i = nA; i < nA + nB; ++i) {
        particles.species[i] = 1;
    }
    const Lattice lattice = cubicLattice(12.0);
    const double volume = 12.0 * 12.0 * 12.0;

    const double productionEnergy = lj.tailEnergyCorrection(particles, lattice);
    const double productionVirial = lj.tailVirialCorrection(particles, lattice);

    // Reference: U_LRC = (2*pi/V) * sum_a sum_b N_a*N_b*I_ab, via the
    // independent quadrature integrals (not the closed-form formula).
    double referenceEnergy = 0.0;
    double referenceVirial = 0.0;
    const std::array<std::pair<std::size_t, LennardJonesParameters>, 2> counted{
        {{nA, species[0]}, {nB, species[1]}}};
    for (const auto& [countA, paramsA] : counted) {
        for (const auto& [countB, paramsB] : counted) {
            const double mixedEpsilon = std::sqrt(paramsA.epsilon * paramsB.epsilon);
            const double mixedSigma = 0.5 * (paramsA.sigma + paramsB.sigma); // Lorentz-Berthelot (default)
            const auto integrals = referenceTailIntegrals(mixedEpsilon, mixedSigma, rc);
            const double weight = static_cast<double>(countA) * static_cast<double>(countB);
            referenceEnergy += weight * integrals.energyIntegral;
            referenceVirial += weight * integrals.virialIntegral;
        }
    }
    referenceEnergy *= 2.0 * std::numbers::pi / volume;
    referenceVirial *= 2.0 * std::numbers::pi / volume;

    REQUIRE(std::abs(productionEnergy - referenceEnergy) / std::abs(referenceEnergy) < 1e-9);
    REQUIRE(std::abs(productionVirial - referenceVirial) / std::abs(referenceVirial) < 1e-9);
}
