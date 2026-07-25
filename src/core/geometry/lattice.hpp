#pragma once

#include <array>

namespace aleator::core {

/// A periodic simulation cell, in the general triclinic (3x3 matrix) form.
///
/// CLAUDE.md invariant #6: triclinic cells are first-class. Many MOFs and
/// zeolites are non-orthogonal, so this type never degrades to three
/// orthogonal lengths + angles internally — `matrix` is always the full
/// lattice-vector matrix, row `i` = lattice vector `a_i` in Å, cartesian
/// components in columns. An orthorhombic cell is simply one where the
/// off-diagonal entries happen to be zero; no code path may assume that.
///
/// All PBC-adjacent operations here are validated against a brute-force
/// ground truth over triclinic and deliberately pathological (e.g.
/// 60/60/60 degree rhombohedral) cells in tests/validation/ — see CLAUDE.md
/// #10 ("Minimum image convention assuming orthorhombic — silently wrong
/// for triclinic MOFs").
class Lattice {
public:
    Lattice() = default;

    /// `matrix[i]` is lattice vector a_i = (matrix[i][0], matrix[i][1], matrix[i][2]), in Å.
    explicit Lattice(const std::array<std::array<double, 3>, 3>& matrix) : matrix_(matrix) {}

    [[nodiscard]] const std::array<std::array<double, 3>, 3>& matrix() const noexcept {
        return matrix_;
    }

    /// Cell volume via the scalar triple product a_0 . (a_1 x a_2), in Å^3.
    /// Pure linear algebra on the stored matrix — not a PBC operation, so
    /// it is implemented directly rather than declared as NotImplemented.
    [[nodiscard]] double volume() const noexcept {
        const auto& m = matrix_;
        const double cx = m[1][1] * m[2][2] - m[1][2] * m[2][1];
        const double cy = m[1][2] * m[2][0] - m[1][0] * m[2][2];
        const double cz = m[1][0] * m[2][1] - m[1][1] * m[2][0];
        return m[0][0] * cx + m[0][1] * cy + m[0][2] * cz;
    }

    /// Converts fractional coordinates (each in, conventionally, [0,1)) to
    /// Cartesian Å: cartesian = fractional * matrix (row-vector convention).
    [[nodiscard]] std::array<double, 3> fractionalToCartesian(
        const std::array<double, 3>& fractional) const;

    /// Converts Cartesian Å coordinates to fractional coordinates, via the
    /// reciprocal-vector formula f_i = cartesian . (a_j x a_k) / volume.
    /// Throws std::invalid_argument if the cell is degenerate (volume ~ 0).
    [[nodiscard]] std::array<double, 3> cartesianToFractional(
        const std::array<double, 3>& cartesian) const;

    /// Minimum-image displacement vector (Å) from point `a` to the nearest
    /// periodic image of `b` under this cell's periodic boundary
    /// conditions. Correct for arbitrary triclinic cells: internally
    /// searches a neighborhood of candidate lattice translations around a
    /// Gauss-reduced basis rather than relying on naive independent
    /// fractional-coordinate rounding, which is wrong for strongly skewed
    /// cells (CLAUDE.md #10). Throws std::invalid_argument if the cell is
    /// degenerate.
    [[nodiscard]] std::array<double, 3> minimumImageDisplacement(
        const std::array<double, 3>& a, const std::array<double, 3>& b) const;

    /// Wraps a Cartesian point into the reference cell (fractional
    /// coordinates in [0,1)^3), returned in Cartesian Å.
    [[nodiscard]] std::array<double, 3> wrapIntoCell(const std::array<double, 3>& cartesian) const;

    /// Perpendicular width (Å) along lattice direction `axis` (0, 1, or 2):
    /// the distance between the two faces of the cell spanned by the other
    /// two lattice vectors, h_axis = volume / |a_j x a_k|. This is the
    /// quantity a minimum-image cutoff must not exceed half of.
    [[nodiscard]] double perpendicularWidth(int axis) const;

    /// Throws std::invalid_argument, naming the offending lattice
    /// direction, if `cutoff` exceeds perpendicularWidth(axis)/2 for any
    /// axis — that configuration lets minimum image double-count the same
    /// neighbor through two different periodic images (CLAUDE.md #10).
    void validateCutoff(double cutoff) const;

    /// Returns an equivalent lattice (same periodicity and volume) whose
    /// basis vectors have been Gauss-reduced: each vector is the shortest
    /// available representative of its coset, which is what makes a small,
    /// bounded neighborhood search sufficient for minimum image on
    /// otherwise badly-conditioned (strongly skewed) cells.
    [[nodiscard]] Lattice reduced() const;

private:
    std::array<std::array<double, 3>, 3> matrix_{};
};

} // namespace aleator::core
