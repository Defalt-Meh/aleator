// Validation anchor for Lattice::minimumImageDisplacement (CLAUDE.md #10:
// "Minimum image convention assuming orthorhombic — silently wrong for
// triclinic MOFs" is a named failure mode). Ground truth here is an
// independent brute-force search over candidate periodic-image
// translations, implemented directly against the raw lattice matrix rather
// than by calling any Lattice method — so a bug shared between the
// production algorithm and the test oracle cannot hide.
//
// Search radius: 27 images (every combination of {-1,0,1} per axis) for
// well-conditioned cells, widened to 125 (every combination of {-2..2}) for
// the deliberately pathological rhombohedral cell, per the milestone spec.

#include <catch2/catch_test_macros.hpp>

#include <array>
#include <cmath>
#include <limits>
#include <utility>

#include "core/geometry/lattice.hpp"

namespace {

using aleator::core::Lattice;
using Vec3 = std::array<double, 3>;

Vec3 sub(const Vec3& a, const Vec3& b) { return {a[0] - b[0], a[1] - b[1], a[2] - b[2]}; }

Vec3 add(const Vec3& a, const Vec3& b) { return {a[0] + b[0], a[1] + b[1], a[2] + b[2]}; }

Vec3 scale(const Vec3& a, double s) { return {a[0] * s, a[1] * s, a[2] * s}; }

double norm(const Vec3& v) { return std::sqrt(v[0] * v[0] + v[1] * v[1] + v[2] * v[2]); }

/// Ground truth, deliberately independent of Lattice's own implementation:
/// enumerate every integer translation n0*a0 + n1*a1 + n2*a2 with each n_i
/// in [-imageRadius, imageRadius], and return the smallest |d - translation|
/// found. `matrix` rows are the raw lattice vectors a0, a1, a2.
double bruteForceMinimumImageDistance(const std::array<std::array<double, 3>, 3>& matrix,
                                       const Vec3& a, const Vec3& b, int imageRadius) {
    const Vec3 a0{matrix[0][0], matrix[0][1], matrix[0][2]};
    const Vec3 a1{matrix[1][0], matrix[1][1], matrix[1][2]};
    const Vec3 a2{matrix[2][0], matrix[2][1], matrix[2][2]};
    const Vec3 d = sub(b, a);

    double best = std::numeric_limits<double>::infinity();
    for (int n0 = -imageRadius; n0 <= imageRadius; ++n0) {
        for (int n1 = -imageRadius; n1 <= imageRadius; ++n1) {
            for (int n2 = -imageRadius; n2 <= imageRadius; ++n2) {
                const Vec3 translation =
                    add(add(scale(a0, static_cast<double>(n0)), scale(a1, static_cast<double>(n1))),
                        scale(a2, static_cast<double>(n2)));
                const double dist = norm(sub(d, translation));
                if (dist < best) {
                    best = dist;
                }
            }
        }
    }
    return best;
}

} // namespace

TEST_CASE("minimumImageDisplacement matches 27-image brute force on a cubic cell",
          "[validation][geometry]") {
    std::array<std::array<double, 3>, 3> m{};
    m[0] = {10.0, 0.0, 0.0};
    m[1] = {0.0, 10.0, 0.0};
    m[2] = {0.0, 0.0, 10.0};
    const Lattice lattice(m);

    const Vec3 a{1.0, 1.0, 1.0};
    const Vec3 b{9.0, 9.0, 9.0}; // nearest image of b is at (-1,-1,-1), distance 2*sqrt(3)

    const auto d = lattice.minimumImageDisplacement(a, b);
    const double got = norm(d);
    const double expected = bruteForceMinimumImageDistance(m, a, b, 1);

    REQUIRE(std::abs(got - expected) < 1e-9);
}

TEST_CASE("minimumImageDisplacement matches 27-image brute force on a moderately triclinic cell",
          "[validation][geometry]") {
    std::array<std::array<double, 3>, 3> m{};
    m[0] = {12.0, 0.0, 0.0};
    m[1] = {3.5, 10.0, 0.0};
    m[2] = {2.0, 1.5, 9.0};
    const Lattice lattice(m);

    // A handful of point pairs, including ones straddling cell boundaries.
    const std::array<std::pair<Vec3, Vec3>, 4> pairs{{
        {Vec3{0.5, 0.5, 0.5}, Vec3{11.0, 9.0, 8.5}},
        {Vec3{1.0, 1.0, 1.0}, Vec3{13.0, 0.5, 0.5}},
        {Vec3{6.0, 5.0, 4.5}, Vec3{6.5, 5.5, 4.0}},
        {Vec3{0.1, 9.9, 8.9}, Vec3{11.9, 0.1, 0.2}},
    }};

    for (const auto& [a, b] : pairs) {
        const auto d = lattice.minimumImageDisplacement(a, b);
        const double got = norm(d);
        const double expected = bruteForceMinimumImageDistance(m, a, b, 1);
        REQUIRE(std::abs(got - expected) < 1e-6);
    }
}

TEST_CASE(
    "minimumImageDisplacement matches 125-image brute force on a pathological 60/60/60 "
    "rhombohedral cell",
    "[validation][geometry]") {
    // A rhombohedral primitive cell with all inter-axis angles at 60
    // degrees (e.g. the FCC primitive cell shape) — exactly the kind of
    // strongly skewed triclinic cell where independently rounding each
    // fractional coordinate (the naive approach) picks the wrong image.
    constexpr double kLen = 8.0;
    constexpr double kCos60 = 0.5;
    std::array<std::array<double, 3>, 3> m{};
    m[0] = {kLen, 0.0, 0.0};
    m[1] = {kLen * kCos60, kLen * std::sqrt(1.0 - kCos60 * kCos60), 0.0};
    // a2 chosen so that a0.a2 = a1.a2 = |a2|^2 * cos(60 deg) = |a2| * kLen * 0.5.
    const double a2z = kLen * std::sqrt(1.0 - 3.0 * kCos60 * kCos60 + 2.0 * kCos60 * kCos60 * kCos60) /
                        std::sqrt(1.0 - kCos60 * kCos60);
    m[2] = {kLen * kCos60, kLen * (kCos60 - kCos60 * kCos60) / std::sqrt(1.0 - kCos60 * kCos60), a2z};
    const Lattice lattice(m);

    // Sanity: this construction really is 60/60/60 between every pair.
    auto dot = [](const Vec3& x, const Vec3& y) { return x[0] * y[0] + x[1] * y[1] + x[2] * y[2]; };
    const Vec3 v0{m[0][0], m[0][1], m[0][2]};
    const Vec3 v1{m[1][0], m[1][1], m[1][2]};
    const Vec3 v2{m[2][0], m[2][1], m[2][2]};
    REQUIRE(std::abs(dot(v0, v1) / (norm(v0) * norm(v1)) - 0.5) < 1e-9);
    REQUIRE(std::abs(dot(v0, v2) / (norm(v0) * norm(v2)) - 0.5) < 1e-9);
    REQUIRE(std::abs(dot(v1, v2) / (norm(v1) * norm(v2)) - 0.5) < 1e-9);

    const std::array<std::pair<Vec3, Vec3>, 3> pairs{{
        {Vec3{0.5, 0.5, 0.5}, Vec3{7.5, 7.0, 6.5}},
        {Vec3{1.0, 0.5, 0.2}, Vec3{6.8, 6.9, 6.7}},
        {Vec3{0.2, 0.2, 0.2}, Vec3{7.8, 7.6, 7.5}},
    }};

    for (const auto& [a, b] : pairs) {
        const auto d = lattice.minimumImageDisplacement(a, b);
        const double got = norm(d);
        const double expected = bruteForceMinimumImageDistance(m, a, b, 2);
        REQUIRE(std::abs(got - expected) < 1e-6);
    }
}
