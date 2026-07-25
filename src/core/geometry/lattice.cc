#include "core/geometry/lattice.hpp"

#include <algorithm>
#include <cmath>
#include <sstream>
#include <stdexcept>

#include "core/exceptions.hpp"

namespace aleator::core {

namespace {

using Vec3 = std::array<double, 3>;

// A volume (or perpendicular area) below this is treated as a degenerate
// cell. Å^3-scale quantities in any physically meaningful cell are many
// orders of magnitude above this, so it only catches genuinely degenerate
// (e.g. default-constructed all-zero) lattices.
constexpr double kDegenerateVolumeThreshold = 1e-10;

Vec3 toVec3(const std::array<double, 3>& row) { return {row[0], row[1], row[2]}; }

Vec3 add(const Vec3& a, const Vec3& b) { return {a[0] + b[0], a[1] + b[1], a[2] + b[2]}; }

Vec3 sub(const Vec3& a, const Vec3& b) { return {a[0] - b[0], a[1] - b[1], a[2] - b[2]}; }

Vec3 scale(const Vec3& a, double s) { return {a[0] * s, a[1] * s, a[2] * s}; }

double dot(const Vec3& a, const Vec3& b) { return a[0] * b[0] + a[1] * b[1] + a[2] * b[2]; }

Vec3 cross(const Vec3& a, const Vec3& b) {
    return {a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0]};
}

double norm2(const Vec3& a) { return dot(a, a); }

std::array<Vec3, 3> rowsOf(const std::array<std::array<double, 3>, 3>& m) {
    return {toVec3(m[0]), toVec3(m[1]), toVec3(m[2])};
}

/// Signed scalar triple product a0 . (a1 x a2) for arbitrary basis vectors
/// (not necessarily this Lattice's own matrix_ — used both for volume() and
/// for the reduced basis produced by greedyReduce()).
double tripleProduct(const Vec3& a0, const Vec3& a1, const Vec3& a2) {
    return dot(a0, cross(a1, a2));
}

/// Fractional coordinates of Cartesian point `c` with respect to basis
/// (b0, b1, b2), via the reciprocal-vector formula f_i = c . (b_j x b_k) /
/// volume. `volume` is passed in rather than recomputed so callers that
/// already have it (e.g. after validating non-degeneracy) don't redo the
/// work.
Vec3 fractionalOf(const Vec3& c, const Vec3& b0, const Vec3& b1, const Vec3& b2, double volume) {
    return {
        dot(c, cross(b1, b2)) / volume,
        dot(c, cross(b2, b0)) / volume,
        dot(c, cross(b0, b1)) / volume,
    };
}

/// Greedy pairwise (Gauss-style) lattice basis reduction: repeatedly
/// shortens each vector by subtracting the nearest integer multiple of
/// another, then re-sorts by length, until a full pass makes no further
/// change (or an iteration cap is hit). This is not a certified Minkowski
/// reducer, but it converges quickly in practice to a basis whose vectors
/// are close to mutually orthogonal in length — which is what lets a
/// bounded neighborhood search (rather than an unbounded one) find the true
/// minimum image even for badly-conditioned (strongly skewed) input cells.
std::array<Vec3, 3> greedyReduce(std::array<Vec3, 3> v) {
    constexpr int kMaxIterations = 50;
    for (int iter = 0; iter < kMaxIterations; ++iter) {
        bool changed = false;
        for (int i = 0; i < 3; ++i) {
            for (int j = 0; j < 3; ++j) {
                if (i == j) {
                    continue;
                }
                const double nj2 = norm2(v[static_cast<std::size_t>(j)]);
                if (nj2 < kDegenerateVolumeThreshold) {
                    continue;
                }
                const double t = std::round(dot(v[static_cast<std::size_t>(i)],
                                                 v[static_cast<std::size_t>(j)]) /
                                             nj2);
                if (t != 0.0) {
                    v[static_cast<std::size_t>(i)] =
                        sub(v[static_cast<std::size_t>(i)],
                            scale(v[static_cast<std::size_t>(j)], t));
                    changed = true;
                }
            }
        }
        std::sort(v.begin(), v.end(),
                  [](const Vec3& x, const Vec3& y) { return norm2(x) < norm2(y); });
        if (!changed) {
            break;
        }
    }
    return v;
}

} // namespace

std::array<double, 3> Lattice::fractionalToCartesian(const std::array<double, 3>& fractional) const {
    const auto& m = matrix_;
    return {
        fractional[0] * m[0][0] + fractional[1] * m[1][0] + fractional[2] * m[2][0],
        fractional[0] * m[0][1] + fractional[1] * m[1][1] + fractional[2] * m[2][1],
        fractional[0] * m[0][2] + fractional[1] * m[1][2] + fractional[2] * m[2][2],
    };
}

std::array<double, 3> Lattice::cartesianToFractional(const std::array<double, 3>& cartesian) const {
    const double vol = volume();
    if (std::abs(vol) < kDegenerateVolumeThreshold) {
        throw std::invalid_argument(
            "Lattice::cartesianToFractional: degenerate lattice (volume ~ 0)");
    }
    const auto [a0, a1, a2] = rowsOf(matrix_);
    const auto f = fractionalOf(toVec3(cartesian), a0, a1, a2, vol);
    return {f[0], f[1], f[2]};
}

std::array<double, 3> Lattice::wrapIntoCell(const std::array<double, 3>& cartesian) const {
    auto frac = cartesianToFractional(cartesian);
    for (double& x : frac) {
        x -= std::floor(x);
    }
    return fractionalToCartesian(frac);
}

std::array<double, 3> Lattice::minimumImageDisplacement(const std::array<double, 3>& a,
                                                          const std::array<double, 3>& b) const {
    const double vol = volume();
    if (std::abs(vol) < kDegenerateVolumeThreshold) {
        throw std::invalid_argument(
            "Lattice::minimumImageDisplacement: degenerate lattice (volume ~ 0)");
    }

    const Vec3 d = sub(toVec3(b), toVec3(a));

    // Search around a Gauss-reduced basis rather than the raw (possibly
    // badly-conditioned) input matrix: reduction is a unimodular change of
    // representative lattice vectors, so it spans exactly the same set of
    // periodic images, but makes a small bounded neighborhood search
    // reliable — naively rounding each fractional coordinate of the
    // *unreduced* matrix independently is the classic wrong shortcut for
    // strongly skewed triclinic cells (CLAUDE.md #10).
    const auto reducedBasis = greedyReduce(rowsOf(matrix_));
    const Vec3& r0 = reducedBasis[0];
    const Vec3& r1 = reducedBasis[1];
    const Vec3& r2 = reducedBasis[2];
    const double reducedVolume = tripleProduct(r0, r1, r2);

    const Vec3 nearest = fractionalOf(d, r0, r1, r2, reducedVolume);
    const auto n0 = std::lround(nearest[0]);
    const auto n1 = std::lround(nearest[1]);
    const auto n2 = std::lround(nearest[2]);

    // A radius-2 (5x5x5 = 125 candidate) neighborhood around the reduced-
    // basis rounding. A certified Minkowski-reduced basis would only need
    // radius 1 (27 candidates); greedyReduce() is a fast heuristic rather
    // than a certified reducer, so the extra margin is a deliberate safety
    // net rather than a proven-necessary cost — verified empirically in
    // tests/validation/test_lattice_minimum_image.cc against a genuine
    // brute-force search, including a deliberately pathological
    // 60/60/60-degree rhombohedral cell.
    constexpr long kSearchRadius = 2;
    Vec3 best = d;
    double bestNorm2 = norm2(d);
    for (long dn0 = -kSearchRadius; dn0 <= kSearchRadius; ++dn0) {
        for (long dn1 = -kSearchRadius; dn1 <= kSearchRadius; ++dn1) {
            for (long dn2 = -kSearchRadius; dn2 <= kSearchRadius; ++dn2) {
                const double c0 = static_cast<double>(n0 + dn0);
                const double c1 = static_cast<double>(n1 + dn1);
                const double c2 = static_cast<double>(n2 + dn2);
                const Vec3 translation =
                    add(add(scale(r0, c0), scale(r1, c1)), scale(r2, c2));
                const Vec3 candidate = sub(d, translation);
                const double candidateNorm2 = norm2(candidate);
                if (candidateNorm2 < bestNorm2) {
                    bestNorm2 = candidateNorm2;
                    best = candidate;
                }
            }
        }
    }
    return {best[0], best[1], best[2]};
}

double Lattice::perpendicularWidth(int axis) const {
    if (axis < 0 || axis > 2) {
        throw std::invalid_argument("Lattice::perpendicularWidth: axis must be 0, 1, or 2");
    }
    const auto [a0, a1, a2] = rowsOf(matrix_);
    Vec3 other0{};
    Vec3 other1{};
    switch (axis) {
        case 0:
            other0 = a1;
            other1 = a2;
            break;
        case 1:
            other0 = a2;
            other1 = a0;
            break;
        default:
            other0 = a0;
            other1 = a1;
            break;
    }
    const double faceArea = std::sqrt(norm2(cross(other0, other1)));
    if (faceArea < kDegenerateVolumeThreshold) {
        throw std::invalid_argument(
            "Lattice::perpendicularWidth: degenerate lattice (zero-area face)");
    }
    return std::abs(volume()) / faceArea;
}

void Lattice::validateCutoff(double cutoff) const {
    static constexpr const char* kAxisNames[3] = {"a0", "a1", "a2"};
    for (int axis = 0; axis < 3; ++axis) {
        const double h = perpendicularWidth(axis);
        if (cutoff > 0.5 * h) {
            std::ostringstream message;
            message << "cutoff (" << cutoff << " Å) exceeds half the perpendicular width "
                    << "along lattice direction " << kAxisNames[axis] << " (L_perp = " << h
                    << " Å, L_perp/2 = " << (0.5 * h) << " Å); minimum image would "
                    << "double-count neighbors through two periodic images -- reduce the cutoff "
                    << "or enlarge the cell along " << kAxisNames[axis] << ".";
            throw std::invalid_argument(message.str());
        }
    }
}

Lattice Lattice::reduced() const {
    const double vol = volume();
    if (std::abs(vol) < kDegenerateVolumeThreshold) {
        throw std::invalid_argument("Lattice::reduced: degenerate lattice (volume ~ 0)");
    }
    const auto r = greedyReduce(rowsOf(matrix_));
    return Lattice({{
        {r[0][0], r[0][1], r[0][2]},
        {r[1][0], r[1][1], r[1][2]},
        {r[2][0], r[2][1], r[2][2]},
    }});
}

} // namespace aleator::core
