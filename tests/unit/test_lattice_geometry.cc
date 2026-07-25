// core/geometry milestone: fractional<->Cartesian round-trip and cutoff
// validation. Minimum-image correctness against a brute-force ground truth
// lives in tests/validation/test_lattice_minimum_image.cc — this file
// covers the more mechanical, unit-level properties.

#include <catch2/catch_test_macros.hpp>

#include <array>
#include <cmath>
#include <cstddef>
#include <stdexcept>
#include <string>

#include "core/geometry/lattice.hpp"

namespace {

using aleator::core::Lattice;

constexpr double kEps = 1e-9;

double norm(const std::array<double, 3>& v) {
    return std::sqrt(v[0] * v[0] + v[1] * v[1] + v[2] * v[2]);
}

std::array<double, 3> sub(const std::array<double, 3>& a, const std::array<double, 3>& b) {
    return {a[0] - b[0], a[1] - b[1], a[2] - b[2]};
}

Lattice cubicLattice(double length) {
    std::array<std::array<double, 3>, 3> m{};
    m[0] = {length, 0.0, 0.0};
    m[1] = {0.0, length, 0.0};
    m[2] = {0.0, 0.0, length};
    return Lattice(m);
}

Lattice triclinicLattice() {
    // A realistic, moderately skewed triclinic cell (not orthogonal, not
    // pathological) — e.g. a typical low-symmetry MOF cell.
    std::array<std::array<double, 3>, 3> m{};
    m[0] = {12.0, 0.0, 0.0};
    m[1] = {3.5, 10.0, 0.0};
    m[2] = {2.0, 1.5, 9.0};
    return Lattice(m);
}

} // namespace

TEST_CASE("fractionalToCartesian and cartesianToFractional round-trip on a cubic cell",
          "[unit][geometry]") {
    const Lattice lattice = cubicLattice(10.0);
    for (const std::array<double, 3>& frac :
         {std::array<double, 3>{0.0, 0.0, 0.0}, std::array<double, 3>{0.25, 0.6, 0.9},
          std::array<double, 3>{-0.3, 1.7, 0.5}}) {
        const auto cart = lattice.fractionalToCartesian(frac);
        const auto back = lattice.cartesianToFractional(cart);
        REQUIRE(norm(sub(back, frac)) < kEps);
    }
}

TEST_CASE("fractionalToCartesian and cartesianToFractional round-trip on a triclinic cell",
          "[unit][geometry]") {
    const Lattice lattice = triclinicLattice();
    for (const std::array<double, 3>& frac :
         {std::array<double, 3>{0.0, 0.0, 0.0}, std::array<double, 3>{0.1, 0.2, 0.3},
          std::array<double, 3>{0.7, 0.4, 0.8}, std::array<double, 3>{-1.2, 2.3, -0.4}}) {
        const auto cart = lattice.fractionalToCartesian(frac);
        const auto back = lattice.cartesianToFractional(cart);
        REQUIRE(norm(sub(back, frac)) < kEps);
    }
}

TEST_CASE("cartesianToFractional throws on a degenerate (zero-volume) lattice",
          "[unit][geometry]") {
    const Lattice degenerate; // default-constructed: all-zero matrix
    REQUIRE_THROWS_AS(degenerate.cartesianToFractional({1.0, 2.0, 3.0}), std::invalid_argument);
}

TEST_CASE("wrapIntoCell brings a Cartesian point back into the reference cell",
          "[unit][geometry]") {
    const Lattice lattice = triclinicLattice();
    const auto cart = lattice.fractionalToCartesian({1.4, -0.6, 2.1});
    const auto wrapped = lattice.wrapIntoCell(cart);
    const auto wrappedFrac = lattice.cartesianToFractional(wrapped);
    for (double f : wrappedFrac) {
        REQUIRE(f >= -kEps);
        REQUIRE(f < 1.0 + kEps);
    }
    // Wrapping must differ from the original point by a whole number of
    // lattice translations, i.e. the fractional parts (mod 1) must match.
    const auto originalFrac = lattice.cartesianToFractional(cart);
    for (std::size_t i = 0; i < 3; ++i) {
        const double diff = wrappedFrac[i] - (originalFrac[i] - std::floor(originalFrac[i]));
        REQUIRE(std::abs(diff) < kEps);
    }
}

TEST_CASE("perpendicularWidth matches h = V/A for a cubic cell", "[unit][geometry]") {
    const Lattice lattice = cubicLattice(10.0);
    for (int axis = 0; axis < 3; ++axis) {
        REQUIRE(std::abs(lattice.perpendicularWidth(axis) - 10.0) < kEps);
    }
}

TEST_CASE("validateCutoff accepts a cutoff within L_perp/2", "[unit][geometry]") {
    const Lattice lattice = cubicLattice(10.0);
    REQUIRE_NOTHROW(lattice.validateCutoff(4.9));
    REQUIRE_NOTHROW(lattice.validateCutoff(5.0));
}

TEST_CASE("validateCutoff rejects a cutoff exceeding L_perp/2 and names the direction",
          "[unit][geometry]") {
    // Elongated along a1 so that direction is the tightest constraint.
    std::array<std::array<double, 3>, 3> m{};
    m[0] = {20.0, 0.0, 0.0};
    m[1] = {0.0, 6.0, 0.0};
    m[2] = {0.0, 0.0, 20.0};
    const Lattice lattice(m);

    REQUIRE_THROWS_AS(lattice.validateCutoff(3.5), std::invalid_argument);
    try {
        lattice.validateCutoff(3.5);
        FAIL("expected validateCutoff to throw");
    } catch (const std::invalid_argument& err) {
        const std::string message = err.what();
        REQUIRE(message.find("a1") != std::string::npos);
    }
}

TEST_CASE("reduced() preserves cell volume", "[unit][geometry]") {
    const Lattice lattice = triclinicLattice();
    const Lattice reduced = lattice.reduced();
    REQUIRE(std::abs(std::abs(reduced.volume()) - std::abs(lattice.volume())) < 1e-6);
}
