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
/// `Lattice` itself only *holds* the cell. Every operation that is a known
/// failure mode if done carelessly (fractional/cartesian conversion, the
/// periodic-boundary minimum-image displacement, wrapping into the
/// reference cell) is declared here as an entry point but not implemented —
/// see CLAUDE.md #10 ("Minimum image convention assuming orthorhombic —
/// silently wrong for triclinic MOFs") and #4 (Ewald / validation anchors).
/// Those are physics-adjacent, correctness-critical algorithms that get
/// their own validation tests (triclinic + cutoff > L/2 rejection) before
/// they are implemented, not stubbed in during scaffolding.
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
    /// Cartesian Å. Declared only — see class-level comment. Throws
    /// NotImplemented.
    [[nodiscard]] std::array<double, 3> fractionalToCartesian(
        const std::array<double, 3>& fractional) const;

    /// Converts Cartesian Å coordinates to fractional coordinates. Declared
    /// only — see class-level comment. Throws NotImplemented.
    [[nodiscard]] std::array<double, 3> cartesianToFractional(
        const std::array<double, 3>& cartesian) const;

    /// Minimum-image displacement vector (Å) from point `a` to point `b`
    /// under this cell's periodic boundary conditions. Must eventually be
    /// correct for triclinic cells and must reject cutoff > L/2 configs at
    /// the caller's input-validation layer, not silently double-count
    /// (CLAUDE.md #10). Declared only. Throws NotImplemented.
    [[nodiscard]] std::array<double, 3> minimumImageDisplacement(
        const std::array<double, 3>& a, const std::array<double, 3>& b) const;

    /// Wraps a Cartesian point into the reference cell [0,1)^3 in fractional
    /// coordinates. Declared only. Throws NotImplemented.
    [[nodiscard]] std::array<double, 3> wrapIntoCell(const std::array<double, 3>& cartesian) const;

private:
    std::array<std::array<double, 3>, 3> matrix_{};
};

} // namespace aleator::core
