#pragma once

#include <array>
#include <string>

namespace aleator::io {

/// An affine symmetry operation on fractional coordinates, as encoded by a
/// CIF `_symmetry_equiv_pos_as_xyz` / `_space_group_symop_operation_xyz`
/// string such as "x,y,z", "-x,-y,-z", or "1/2+x,1/2-y,z":
///
///     x' = matrix[0][0]*x + matrix[0][1]*y + matrix[0][2]*z + translation[0]
///     y' = matrix[1][0]*x + matrix[1][1]*y + matrix[1][2]*z + translation[1]
///     z' = matrix[2][0]*x + matrix[2][1]*y + matrix[2][2]*z + translation[2]
struct SymmetryOperation {
    std::array<std::array<double, 3>, 3> matrix{};
    std::array<double, 3> translation{};

    [[nodiscard]] std::array<double, 3> apply(const std::array<double, 3>& fractional) const;
};

/// Parses one comma-separated CIF symmetry-equivalent-position string
/// ("x,y,z", "-y,x-y,z+1/3", ...) into a SymmetryOperation. Coefficients on
/// x/y/z and constant terms may be integers, decimals, or simple fractions
/// (e.g. "1/2", "2/3"). Throws CifParseError on malformed input (wrong
/// number of components, or a component that isn't a sum of signed
/// x/y/z terms and/or a signed constant).
[[nodiscard]] SymmetryOperation parseSymmetryOperation(const std::string& text);

} // namespace aleator::io
