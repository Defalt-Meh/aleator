// Validation anchor (CLAUDE.md #4: "MD forces | Analytic force vs. central
// finite difference of energy | 1e-6"). Exercises the full pipeline —
// VerletList, minimum image, species mixing — not just the isolated
// pairwise formula: for each particle and each coordinate, perturb by +-h
// and compare -dE/dcoord (central difference) to the analytic force
// component from computeForces(), for every truncation scheme.

#include <catch2/catch_test_macros.hpp>

#include <cmath>
#include <random>
#include <vector>

#include "core/geometry/lattice.hpp"
#include "core/math/particle_data.hpp"
#include "core/memory/aligned_allocator.hpp"
#include "core/neighbor/verlet_list.hpp"
#include "forcefield/pairwise/lennard_jones.hpp"

namespace {

using aleator::core::Lattice;
using aleator::core::ParticleData;
using aleator::core::VerletList;
using aleator::forcefield::Forces;
using aleator::forcefield::LennardJones;
using aleator::forcefield::LennardJonesParameters;
using aleator::forcefield::LennardJonesTruncation;
using aleator::forcefield::MixingRule;

Lattice cubicLattice(double length) {
    std::array<std::array<double, 3>, 3> m{};
    m[0] = {length, 0.0, 0.0};
    m[1] = {0.0, length, 0.0};
    m[2] = {0.0, 0.0, length};
    return Lattice(m);
}

/// Places particles on a simple-cubic lattice with small random jitter,
/// rather than uniformly at random: a central finite difference is only
/// accurate where the energy landscape is smooth on the scale of the step
/// size, and the LJ repulsive wall (~1/r^12) is emphatically not smooth
/// near r -> 0. Fully random placement occasionally puts two particles
/// pathologically close together by chance, which blows up the finite
/// difference's truncation error (a real numerical artifact, not a force
/// bug) — a lattice + small jitter keeps every pair at a sane separation
/// while still exercising a range of distances and species combinations.
ParticleData randomParticles(std::size_t count, double boxLength, unsigned speciesCount,
                              unsigned seed) {
    std::mt19937 rng(seed);
    std::uniform_int_distribution<unsigned> speciesDist(0, speciesCount - 1);

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
                particles.species[placed] = speciesDist(rng);
                ++placed;
            }
        }
    }
    return particles;
}

/// Checks -dE/dcoord (central finite difference, step h) against the
/// analytic force from computeForces(), for every particle/axis.
void checkForcesAgainstFiniteDifference(const LennardJones& lj, ParticleData particles,
                                         const Lattice& lattice, double cutoff, double skin,
                                         double h, double tolerance) {
    VerletList neighbors(skin);
    neighbors.build(particles, lattice, cutoff);

    Forces analytic;
    analytic.resize(particles.size());
    lj.computeForces(particles, lattice, neighbors, analytic);

    auto perturbedEnergy = [&](std::vector<double, aleator::core::AlignedAllocator<double>>& coordArray,
                                std::size_t i, double delta) {
        const double original = coordArray[i];
        coordArray[i] = original + delta;
        const double energy = lj.computeEnergy(particles, lattice, neighbors);
        coordArray[i] = original;
        return energy;
    };

    for (std::size_t i = 0; i < particles.size(); ++i) {
        const double dEdx = (perturbedEnergy(particles.x, i, h) - perturbedEnergy(particles.x, i, -h)) /
                             (2.0 * h);
        const double dEdy = (perturbedEnergy(particles.y, i, h) - perturbedEnergy(particles.y, i, -h)) /
                             (2.0 * h);
        const double dEdz = (perturbedEnergy(particles.z, i, h) - perturbedEnergy(particles.z, i, -h)) /
                             (2.0 * h);

        INFO("particle " << i);
        REQUIRE(std::abs(analytic.fx[i] - (-dEdx)) < tolerance);
        REQUIRE(std::abs(analytic.fy[i] - (-dEdy)) < tolerance);
        REQUIRE(std::abs(analytic.fz[i] - (-dEdz)) < tolerance);
    }
}

} // namespace

TEST_CASE("LennardJones analytic forces match central finite difference: Truncated -- "
          "single species",
          "[validation][forcefield]") {
    const Lattice lattice = cubicLattice(10.0);
    const ParticleData particles = randomParticles(40, 10.0, 1, 1);
    std::vector<LennardJonesParameters> species{{1.2, 1.0}};
    LennardJones lj(species, 3.0, LennardJonesTruncation::Truncated);
    checkForcesAgainstFiniteDifference(lj, particles, lattice, 3.0, /*skin=*/0.8, /*h=*/1e-6,
                                        /*tolerance=*/1e-6);
}

TEST_CASE("LennardJones analytic forces match central finite difference: Shifted -- "
          "single species",
          "[validation][forcefield]") {
    const Lattice lattice = cubicLattice(10.0);
    const ParticleData particles = randomParticles(40, 10.0, 1, 2);
    std::vector<LennardJonesParameters> species{{1.0, 1.0}};
    LennardJones lj(species, 3.0, LennardJonesTruncation::Shifted);
    checkForcesAgainstFiniteDifference(lj, particles, lattice, 3.0, /*skin=*/0.8, /*h=*/1e-6,
                                        /*tolerance=*/1e-6);
}

TEST_CASE("LennardJones analytic forces match central finite difference: LinearForceShifted -- "
          "single species",
          "[validation][forcefield]") {
    const Lattice lattice = cubicLattice(10.0);
    const ParticleData particles = randomParticles(40, 10.0, 1, 3);
    std::vector<LennardJonesParameters> species{{1.0, 1.0}};
    LennardJones lj(species, 3.0, LennardJonesTruncation::LinearForceShifted);
    checkForcesAgainstFiniteDifference(lj, particles, lattice, 3.0, /*skin=*/0.8, /*h=*/1e-6,
                                        /*tolerance=*/1e-6);
}

TEST_CASE("LennardJones analytic forces match central finite difference: two species -- "
          "Lorentz-Berthelot mixing",
          "[validation][forcefield]") {
    const Lattice lattice = cubicLattice(10.0);
    const ParticleData particles = randomParticles(50, 10.0, 2, 4);
    std::vector<LennardJonesParameters> species{{1.0, 1.0}, {0.6, 1.4}};
    LennardJones lj(species, 3.0, LennardJonesTruncation::LinearForceShifted,
                     MixingRule::LorentzBerthelot);
    checkForcesAgainstFiniteDifference(lj, particles, lattice, 3.0, /*skin=*/0.8, /*h=*/1e-6,
                                        /*tolerance=*/1e-6);
}

TEST_CASE("LennardJones analytic forces match central finite difference: two species -- "
          "geometric mixing",
          "[validation][forcefield]") {
    const Lattice lattice = cubicLattice(10.0);
    const ParticleData particles = randomParticles(50, 10.0, 2, 5);
    std::vector<LennardJonesParameters> species{{1.0, 1.0}, {2.0, 0.8}};
    LennardJones lj(species, 3.0, LennardJonesTruncation::Truncated, MixingRule::Geometric);
    checkForcesAgainstFiniteDifference(lj, particles, lattice, 3.0, /*skin=*/0.8, /*h=*/1e-6,
                                        /*tolerance=*/1e-6);
}
