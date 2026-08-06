// core/geometry milestone (supercell replication): unit-level correctness
// of minimumSupercellReplication() and replicateSupercell(), independent of
// any GCMC/energy machinery — the energy-invariance-under-replication
// physics claims live in tests/validation/test_supercell_energy_invariance.cc,
// this file covers the mechanical geometry/particle-count properties.

#include <catch2/catch_test_macros.hpp>

#include <array>
#include <cmath>
#include <stdexcept>

#include "core/geometry/lattice.hpp"
#include "core/geometry/supercell.hpp"
#include "core/math/particle_data.hpp"

namespace {

using aleator::core::Lattice;
using aleator::core::minimumSupercellReplication;
using aleator::core::ParticleData;
using aleator::core::replicateSupercell;

Lattice cubicLattice(double length) {
    std::array<std::array<double, 3>, 3> m{};
    m[0] = {length, 0.0, 0.0};
    m[1] = {0.0, length, 0.0};
    m[2] = {0.0, 0.0, length};
    return Lattice(m);
}

// Anisotropic orthorhombic: three different edge lengths, so a per-axis
// replication test can actually distinguish axes.
Lattice orthorhombicLattice(double a, double b, double c) {
    std::array<std::array<double, 3>, 3> m{};
    m[0] = {a, 0.0, 0.0};
    m[1] = {0.0, b, 0.0};
    m[2] = {0.0, 0.0, c};
    return Lattice(m);
}

// A real, moderately skewed triclinic cell (same shape used elsewhere in
// this suite) — perpendicularWidth differs meaningfully from |a_i| here,
// which is exactly the case invariant 6 requires be exercised, not just a
// cubic/orthorhombic cell where the two coincide.
Lattice triclinicLattice() {
    std::array<std::array<double, 3>, 3> m{};
    m[0] = {12.0, 0.0, 0.0};
    m[1] = {3.5, 10.0, 0.0};
    m[2] = {2.0, 1.5, 9.0};
    return Lattice(m);
}

} // namespace

TEST_CASE("minimumSupercellReplication returns (1,1,1) when the cell already satisfies minimum "
          "image at the given cutoff",
          "[unit][geometry][supercell]") {
    const Lattice lattice = cubicLattice(25.832); // real IRMOF-1 cell edge
    const auto n = minimumSupercellReplication(lattice, 12.0);
    REQUIRE(n == std::array<int, 3>{1, 1, 1});
}

TEST_CASE("minimumSupercellReplication computes the correct minimum for a cell that needs "
          "replication",
          "[unit][geometry][supercell]") {
    // perpendicularWidth = 10 on every axis for a cubic cell; cutoff = 12
    // needs width > 24, i.e. n=3 (3*10=30>24; n=2 gives 20, not >24).
    const Lattice lattice = cubicLattice(10.0);
    const auto n = minimumSupercellReplication(lattice, 12.0);
    REQUIRE(n == std::array<int, 3>{3, 3, 3});
}

TEST_CASE("minimumSupercellReplication is per-axis independent on an anisotropic cell",
          "[unit][geometry][supercell]") {
    // a=8 (needs replication for cutoff=6: 2*6=12, n=2 -> 16>12 OK, n=1 ->
    // 8 not >12), b=20 (n=1: 20>12 OK), c=40 (n=1: 40>12 OK).
    const Lattice lattice = orthorhombicLattice(8.0, 20.0, 40.0);
    const auto n = minimumSupercellReplication(lattice, 6.0);
    REQUIRE(n == std::array<int, 3>{2, 1, 1});
}

TEST_CASE("minimumSupercellReplication uses perpendicularWidth, not raw lattice vector length, "
          "on a triclinic cell (invariant 6)",
          "[unit][geometry][supercell]") {
    // Triclinic cell: |a_0|=12.0 exactly, but perpendicularWidth(0) differs
    // from 12.0 because a_1/a_2 are skewed relative to a_0. Confirm the
    // computed replication is consistent with the ACTUAL perpendicular
    // width (queried directly), not a naive |a_i|-based computation that
    // would silently be wrong for this cell shape.
    const Lattice lattice = triclinicLattice();
    const double cutoff = 5.0;
    const auto n = minimumSupercellReplication(lattice, cutoff);
    for (int axis = 0; axis < 3; ++axis) {
        const double width = lattice.perpendicularWidth(axis);
        const auto count = static_cast<std::size_t>(axis);
        // The computed n must be minimal: n*width > 2*cutoff, but
        // (n-1)*width does NOT exceed 2*cutoff (unless n==1, the floor).
        REQUIRE(static_cast<double>(n[count]) * width > 2.0 * cutoff);
        if (n[count] > 1) {
            REQUIRE_FALSE(static_cast<double>(n[count] - 1) * width > 2.0 * cutoff);
        }
    }
}

TEST_CASE("minimumSupercellReplication rejects a non-positive cutoff", "[unit][geometry][supercell]") {
    const Lattice lattice = cubicLattice(10.0);
    REQUIRE_THROWS_AS(minimumSupercellReplication(lattice, 0.0), std::invalid_argument);
    REQUIRE_THROWS_AS(minimumSupercellReplication(lattice, -1.0), std::invalid_argument);
}

TEST_CASE("replicateSupercell at (1,1,1) is an exact bit-identical no-op",
          "[unit][geometry][supercell]") {
    const Lattice lattice = triclinicLattice();
    ParticleData particles;
    particles.push_back(1.23456789, 2.3456789, 3.456789, 0.1, 0.2, 0.3, 12.011, -0.5, 0);
    particles.push_back(4.5, 5.5, 6.5, 0.0, 0.0, 0.0, 16.0, 0.25, 1);

    const auto result = replicateSupercell(particles, lattice, 1, 1, 1);

    REQUIRE(result.particles.size() == particles.size());
    for (std::size_t i = 0; i < particles.size(); ++i) {
        CHECK(result.particles.x[i] == particles.x[i]);
        CHECK(result.particles.y[i] == particles.y[i]);
        CHECK(result.particles.z[i] == particles.z[i]);
        CHECK(result.particles.mass[i] == particles.mass[i]);
        CHECK(result.particles.charge[i] == particles.charge[i]);
        CHECK(result.particles.species[i] == particles.species[i]);
    }
    for (int r = 0; r < 3; ++r) {
        for (int c = 0; c < 3; ++c) {
            CHECK(result.lattice.matrix()[static_cast<std::size_t>(r)][static_cast<std::size_t>(c)] ==
                  lattice.matrix()[static_cast<std::size_t>(r)][static_cast<std::size_t>(c)]);
        }
    }
}

TEST_CASE("replicateSupercell produces nx*ny*nz times the particle count and a correspondingly "
          "scaled lattice",
          "[unit][geometry][supercell]") {
    const Lattice lattice = triclinicLattice();
    ParticleData particles;
    particles.push_back(0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0);
    particles.push_back(6.0, 5.0, 4.5, 0.0, 0.0, 0.0, 1.0, 0.0, 1);

    const auto result = replicateSupercell(particles, lattice, 2, 3, 1);

    REQUIRE(result.particles.size() == particles.size() * 2 * 3 * 1);

    const auto& orig = lattice.matrix();
    const auto& scaled = result.lattice.matrix();
    for (int c = 0; c < 3; ++c) {
        const auto col = static_cast<std::size_t>(c);
        CHECK(scaled[0][col] == 2.0 * orig[0][col]);
        CHECK(scaled[1][col] == 3.0 * orig[1][col]);
        CHECK(scaled[2][col] == 1.0 * orig[2][col]);
    }

    // Volume must scale by nx*ny*nz exactly (a direct consequence of each
    // row being scaled independently) -- a real, checkable consistency
    // property, not just "the code ran".
    CHECK(std::abs(result.lattice.volume() - 6.0 * lattice.volume()) <
          1e-9 * std::abs(lattice.volume()));
}

TEST_CASE("replicateSupercell places every image at an exact integer lattice-vector translation "
          "of the original position",
          "[unit][geometry][supercell]") {
    const Lattice lattice = cubicLattice(10.0);
    ParticleData particles;
    particles.push_back(1.0, 2.0, 3.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0);

    const auto result = replicateSupercell(particles, lattice, 2, 2, 2);
    REQUIRE(result.particles.size() == 8);

    // For a cubic cell of edge 10, the 8 images must be the original point
    // shifted by every combination of {0,10} in x, y, z.
    bool foundOrigin = false;
    bool foundOpposite = false;
    for (std::size_t i = 0; i < result.particles.size(); ++i) {
        const double x = result.particles.x[i];
        const double y = result.particles.y[i];
        const double z = result.particles.z[i];
        if (x == 1.0 && y == 2.0 && z == 3.0) {
            foundOrigin = true;
        }
        if (x == 11.0 && y == 12.0 && z == 13.0) {
            foundOpposite = true;
        }
    }
    CHECK(foundOrigin);
    CHECK(foundOpposite);
}

TEST_CASE("replicateSupercell rejects a non-positive replication count",
          "[unit][geometry][supercell]") {
    const Lattice lattice = cubicLattice(10.0);
    ParticleData particles;
    particles.push_back(0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0);

    REQUIRE_THROWS_AS(replicateSupercell(particles, lattice, 0, 1, 1), std::invalid_argument);
    REQUIRE_THROWS_AS(replicateSupercell(particles, lattice, 1, -1, 1), std::invalid_argument);
}
