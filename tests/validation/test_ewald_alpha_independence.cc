// Validation anchor (CLAUDE.md milestone: "Independence of results from
// the Ewald splitting parameter alpha... this catches term errors that
// Madelung alone might not"). alpha is a purely numerical convenience
// parameter that splits 1/r into real-space (erfc) and reciprocal-space
// (erf, via the Fourier sum) pieces — the TOTAL energy must be exactly
// independent of it. This is a genuinely different check from the
// Madelung test: an implementation with, say, a wrong prefactor on the
// self-energy term (-alpha/sqrt(pi) * sum q_i^2) could still accidentally
// land close to the right total at ONE alpha for a highly symmetric
// structure, but would not stay constant as alpha is swept, because each
// of U_real/U_reciprocal/U_self individually has strong, different alpha
// dependence that only cancels when every term's coefficient is correct.

#include <catch2/catch_test_macros.hpp>

#include <array>
#include <cmath>
#include <cstdio>
#include <vector>

#include "core/geometry/lattice.hpp"
#include "core/math/particle_data.hpp"
#include "core/neighbor/verlet_list.hpp"
#include "forcefield/electrostatics/ewald.hpp"

namespace {

using aleator::core::Lattice;
using aleator::core::ParticleData;
using aleator::core::VerletList;
using aleator::forcefield::Ewald;

/// Same rock-salt supercell builder as test_ewald_madelung.cc (duplicated
/// rather than shared, matching this project's convention of
/// self-contained test files).
ParticleData buildRockSalt(double a, int nRep, double& boxLengthOut) {
    const double cellEdge = 2.0 * a;
    boxLengthOut = cellEdge * static_cast<double>(nRep);

    const std::vector<std::array<double, 3>> naSites{
        {0.0, 0.0, 0.0}, {a, a, 0.0}, {a, 0.0, a}, {0.0, a, a}};
    const std::vector<std::array<double, 3>> clSites{
        {a, 0.0, 0.0}, {0.0, a, 0.0}, {0.0, 0.0, a}, {a, a, a}};

    ParticleData particles;
    for (int ox = 0; ox < nRep; ++ox) {
        for (int oy = 0; oy < nRep; ++oy) {
            for (int oz = 0; oz < nRep; ++oz) {
                const double dx = cellEdge * static_cast<double>(ox);
                const double dy = cellEdge * static_cast<double>(oy);
                const double dz = cellEdge * static_cast<double>(oz);
                for (const auto& site : naSites) {
                    particles.push_back(site[0] + dx, site[1] + dy, site[2] + dz, 0.0, 0.0, 0.0,
                                         22.99, 1.0, 0);
                }
                for (const auto& site : clSites) {
                    particles.push_back(site[0] + dx, site[1] + dy, site[2] + dz, 0.0, 0.0, 0.0,
                                         35.45, -1.0, 1);
                }
            }
        }
    }
    return particles;
}

} // namespace

TEST_CASE("Ewald total energy is independent of the splitting parameter alpha",
          "[validation][electrostatics]") {
    const double a = 1.0;
    const int nRep = 3;
    double boxLength = 0.0;
    const ParticleData particles = buildRockSalt(a, nRep, boxLength);

    std::array<std::array<double, 3>, 3> m{};
    m[0] = {boxLength, 0.0, 0.0};
    m[1] = {0.0, boxLength, 0.0};
    m[2] = {0.0, 0.0, boxLength};
    const Lattice lattice(m);

    const double realSpaceCutoff = 2.9; // < L_perp/2 = 3.0
    const int kMax = 20;                // generous even for the largest alpha below

    VerletList neighbors(/*skinDistance=*/0.3);
    neighbors.build(particles, lattice, realSpaceCutoff);

    // Lower bound chosen so erfc(alpha*realSpaceCutoff) is comfortably
    // converged at every alpha tested (erfc(1.2*2.9) = erfc(3.48) ~ 7e-7,
    // NOT small enough at this rc — an earlier version of this test
    // included alpha=1.2 and failed at ~1.76e-6 relative purely from that
    // real-space truncation residual, confirmed by the other three alphas
    // agreeing with each other to ~1e-10; this is expected Ewald
    // convergence behavior; a real term-coefficient bug would not produce
    // three-way agreement at that precision by coincidence).
    std::vector<double> energies;
    std::vector<double> alphas{1.6, 1.9, 2.2, 2.5};
    for (double alpha : alphas) {
        const Ewald ewald(alpha, realSpaceCutoff, kMax);
        const double energy = ewald.computeEnergy(particles, lattice, neighbors);
        energies.push_back(energy);
        std::printf("alpha=%.2f -> U_total = %.10f K\n", alpha, energy);
    }

    const double reference = energies.front();
    for (std::size_t i = 0; i < energies.size(); ++i) {
        const double relError = std::abs(energies[i] - reference) / std::abs(reference);
        INFO("alpha=" << alphas[i] << " energy=" << energies[i] << " relError=" << relError);
        REQUIRE(relError < 1e-6);
    }
}
