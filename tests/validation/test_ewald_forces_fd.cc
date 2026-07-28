// Validation anchor (CLAUDE.md #4: "MD forces | Analytic force vs. central
// finite difference of energy | 1e-6"), applied to Ewald specifically:
// checks the real-space, reciprocal-space, AND exclusion-correction force
// terms together, since all three contribute to the total energy's
// gradient.

#include <catch2/catch_test_macros.hpp>

#include <algorithm>
#include <array>
#include <cmath>
#include <random>
#include <vector>

#include "core/geometry/lattice.hpp"
#include "core/math/particle_data.hpp"
#include "core/memory/aligned_allocator.hpp"
#include "core/neighbor/verlet_list.hpp"
#include "forcefield/electrostatics/ewald.hpp"

namespace {

using aleator::core::Lattice;
using aleator::core::ParticleData;
using aleator::core::VerletList;
using aleator::forcefield::Ewald;
using aleator::forcefield::Forces;

Lattice cubicLattice(double length) {
    std::array<std::array<double, 3>, 3> m{};
    m[0] = {length, 0.0, 0.0};
    m[1] = {0.0, length, 0.0};
    m[2] = {0.0, 0.0, length};
    return Lattice(m);
}

/// Jittered-lattice placement (see the analogous helper and comment in
/// tests/validation/test_lennard_jones_forces_fd.cc — fully random
/// placement risks pathologically close pairs that blow up finite-
/// difference truncation error, unrelated to any force bug). Alternates
/// +1/-1 charges so the system is exactly neutral.
ParticleData neutralJitteredParticles(std::size_t count, double boxLength, unsigned seed) {
    std::mt19937 rng(seed);
    const auto gridSize = static_cast<std::size_t>(std::ceil(std::cbrt(static_cast<double>(count))));
    const double spacing = boxLength / static_cast<double>(gridSize);
    std::uniform_real_distribution<double> jitter(-0.15 * spacing, 0.15 * spacing);

    ParticleData particles;
    particles.resize(count);
    std::size_t placed = 0;
    for (std::size_t ix = 0; ix < gridSize && placed < count; ++ix) {
        for (std::size_t iy = 0; iy < gridSize && placed < count; ++iy) {
            for (std::size_t iz = 0; iz < gridSize && placed < count; ++iz) {
                particles.x[placed] = (static_cast<double>(ix) + 0.5) * spacing + jitter(rng);
                particles.y[placed] = (static_cast<double>(iy) + 0.5) * spacing + jitter(rng);
                particles.z[placed] = (static_cast<double>(iz) + 0.5) * spacing + jitter(rng);
                particles.charge[placed] = (placed % 2 == 0) ? 1.0 : -1.0;
                particles.species[placed] = static_cast<std::uint32_t>(placed % 2);
                ++placed;
            }
        }
    }
    return particles;
}

/// Ewald forces live on a much larger absolute scale than the LJ
/// reduced-unit forces checked elsewhere (Ewald::kCoulombConstant ~
/// 1.67e5 K*Angstrom, so typical force components here are O(1e4-1e5)
/// K/Angstrom) — a flat 1e-6 *absolute* tolerance would really be
/// demanding ~1e-11 *relative* precision, tighter than a central
/// difference at h=1e-6 can deliver given ordinary floating-point
/// roundoff at that magnitude. CLAUDE.md's "1e-6" MD-forces tolerance is
/// applied here the same way the "1e-8 relative" LJ-energy tolerance is:
/// as a relative figure, matching the spirit of the other validation
/// anchors in that table. A small absolute floor (in K/Angstrom) handles
/// force components that are near zero, where a relative check alone
/// would be meaningless.
bool closeEnough(double analytic, double finiteDifference, double relativeTolerance = 1e-6,
                  double absoluteFloor = 1e-3) {
    const double diff = std::abs(analytic - finiteDifference);
    if (diff < absoluteFloor) {
        return true;
    }
    return diff / std::max(std::abs(analytic), std::abs(finiteDifference)) < relativeTolerance;
}

} // namespace

TEST_CASE("Ewald analytic forces match central finite difference of the energy",
          "[validation][electrostatics]") {
    const double boxLength = 6.0;
    const Lattice lattice = cubicLattice(boxLength);
    ParticleData particles = neutralJitteredParticles(16, boxLength, 7);

    const double alpha = 1.8;
    const double realSpaceCutoff = 2.8; // < L_perp/2 = 3.0
    const int kMax = 16;

    VerletList neighbors(/*skinDistance=*/0.8);
    neighbors.build(particles, lattice, realSpaceCutoff);

    const Ewald ewald(alpha, realSpaceCutoff, kMax);

    Forces analytic;
    analytic.resize(particles.size());
    ewald.computeForces(particles, lattice, neighbors, analytic);

    const double h = 1e-4;

    auto perturbedEnergy =
        [&](std::vector<double, aleator::core::AlignedAllocator<double>>& coordArray, std::size_t i,
            double delta) {
            const double original = coordArray[i];
            coordArray[i] = original + delta;
            const double energy = ewald.computeEnergy(particles, lattice, neighbors);
            coordArray[i] = original;
            return energy;
        };

    for (std::size_t i = 0; i < particles.size(); ++i) {
        const double dEdx =
            (perturbedEnergy(particles.x, i, h) - perturbedEnergy(particles.x, i, -h)) /
            (2.0 * h);
        const double dEdy =
            (perturbedEnergy(particles.y, i, h) - perturbedEnergy(particles.y, i, -h)) /
            (2.0 * h);
        const double dEdz =
            (perturbedEnergy(particles.z, i, h) - perturbedEnergy(particles.z, i, -h)) /
            (2.0 * h);

        INFO("particle " << i << " analytic=(" << analytic.fx[i] << "," << analytic.fy[i] << ","
                          << analytic.fz[i] << ") fd=(" << -dEdx << "," << -dEdy << "," << -dEdz
                          << ")");
        CHECK(closeEnough(analytic.fx[i], -dEdx));
        CHECK(closeEnough(analytic.fy[i], -dEdy));
        CHECK(closeEnough(analytic.fz[i], -dEdz));
    }
}

TEST_CASE("Ewald analytic forces match central finite difference with intramolecular exclusions",
          "[validation][electrostatics]") {
    // 4 "molecules" of 2 same-sign-adjacent... no: pairs of oppositely
    // charged atoms (i, i+1) bonded (excluded from direct Coulomb), each
    // molecule itself neutral, exercising the exclusion-correction force
    // term specifically.
    const double boxLength = 6.0;
    const Lattice lattice = cubicLattice(boxLength);
    ParticleData particles = neutralJitteredParticles(16, boxLength, 11);

    std::vector<std::pair<std::size_t, std::size_t>> exclusions;
    for (std::size_t i = 0; i + 1 < particles.size(); i += 2) {
        exclusions.emplace_back(i, i + 1);
    }

    const double alpha = 1.8;
    const double realSpaceCutoff = 2.8;
    const int kMax = 16;

    VerletList neighbors(/*skinDistance=*/0.8);
    neighbors.build(particles, lattice, realSpaceCutoff);

    const Ewald ewald(alpha, realSpaceCutoff, kMax, exclusions);

    Forces analytic;
    analytic.resize(particles.size());
    ewald.computeForces(particles, lattice, neighbors, analytic);

    const double h = 1e-4;
    auto perturbedEnergy =
        [&](std::vector<double, aleator::core::AlignedAllocator<double>>& coordArray, std::size_t i,
            double delta) {
            const double original = coordArray[i];
            coordArray[i] = original + delta;
            const double energy = ewald.computeEnergy(particles, lattice, neighbors);
            coordArray[i] = original;
            return energy;
        };

    for (std::size_t i = 0; i < particles.size(); ++i) {
        const double dEdx =
            (perturbedEnergy(particles.x, i, h) - perturbedEnergy(particles.x, i, -h)) / (2.0 * h);
        const double dEdy =
            (perturbedEnergy(particles.y, i, h) - perturbedEnergy(particles.y, i, -h)) / (2.0 * h);
        const double dEdz =
            (perturbedEnergy(particles.z, i, h) - perturbedEnergy(particles.z, i, -h)) / (2.0 * h);

        INFO("particle " << i << " analytic=(" << analytic.fx[i] << "," << analytic.fy[i] << ","
                          << analytic.fz[i] << ") fd=(" << -dEdx << "," << -dEdy << "," << -dEdz
                          << ")");
        CHECK(closeEnough(analytic.fx[i], -dEdx));
        CHECK(closeEnough(analytic.fy[i], -dEdy));
        CHECK(closeEnough(analytic.fz[i], -dEdz));
    }
}
