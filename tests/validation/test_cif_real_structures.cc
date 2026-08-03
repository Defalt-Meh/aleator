// Validation anchor for the CIF reader (CLAUDE.md io/ milestone): parses
// real structure files from the IZA Database of Zeolite Structures
// (Baerlocher & McCusker, http://www.iza-structure.org/databases/,
// fetched verbatim 2026-07 as LTA.cif and PTY.cif) and checks the result
// against independently-derived ground truth, not just "did it not crash".
//
// LTA: cubic, space group Pm-3m (#221) — 48 symmetry operations, the
// "non-trivial symmetry" case. PTY: triclinic, space group P-1 (#2) — the
// triclinic case, and small enough (2 symmetry operations, 16
// asymmetric-unit atoms) to hand-verify every expected multiplicity by
// site-symmetry analysis rather than trusting a remembered database number
// (see the comment above that test for the derivation).
//
// Every zeolite framework is a 4-connected net of corner-sharing TO4
// tetrahedra, so #O == 2 * #T always holds after correct symmetry
// expansion regardless of space group — that invariant is checked for both
// structures as a space-group-independent correctness proof, in addition
// to the specific atom counts.

#include <catch2/catch_test_macros.hpp>

#include <array>
#include <cmath>
#include <cstdio>
#include <filesystem>
#include <fstream>
#include <numbers>
#include <string>

#include "io/cif_error.hpp"
#include "io/structure_io.hpp"

namespace {

using aleator::io::readCif;
using aleator::io::StructureData;

/// Writes `content` to a temp file and returns its path; the file is
/// removed when the returned guard goes out of scope.
class TempCifFile {
public:
    TempCifFile(const std::string& name, const std::string& content)
        : path_(std::filesystem::temp_directory_path() / name) {
        std::ofstream out(path_);
        out << content;
    }
    ~TempCifFile() { std::filesystem::remove(path_); }
    TempCifFile(const TempCifFile&) = delete;
    TempCifFile& operator=(const TempCifFile&) = delete;

    [[nodiscard]] const std::filesystem::path& path() const { return path_; }

private:
    std::filesystem::path path_;
};

double vectorLength(const std::array<double, 3>& v) {
    return std::sqrt(v[0] * v[0] + v[1] * v[1] + v[2] * v[2]);
}

double angleBetweenDegrees(const std::array<double, 3>& a, const std::array<double, 3>& b) {
    const double dot = a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
    const double cosine = dot / (vectorLength(a) * vectorLength(b));
    return std::acos(cosine) * 180.0 / std::numbers::pi;
}

int countSpecies(const StructureData& structure, const std::string& element) {
    int index = -1;
    for (std::size_t i = 0; i < structure.speciesSymbols.size(); ++i) {
        if (structure.speciesSymbols[i] == element) {
            index = static_cast<int>(i);
            break;
        }
    }
    if (index < 0) {
        return 0;
    }
    int count = 0;
    for (auto species : structure.particles.species) {
        if (static_cast<int>(species) == index) {
            ++count;
        }
    }
    return count;
}

// clang-format off
constexpr const char* kLtaCif = R"CIF(data_LTA
#**************************************************************************
#
# CIF taken from the IZA-SC Database of Zeolite Structures
# Ch. Baerlocher and L.B. McCusker
# Database of Zeolite Structures: http://www.iza-structure.org/databases/
#
# The atom coordinates and the cell parameters were optimized with DLS76
# assuming a pure SiO2 composition.
#
#**************************************************************************

_cell_length_a                  11.9190(0)
_cell_length_b                  11.9190(0)
_cell_length_c                  11.9190(0)
_cell_angle_alpha               90.0000(0)
_cell_angle_beta                90.0000(0)
_cell_angle_gamma               90.0000(0)

_symmetry_space_group_name_H-M     'P m -3 m'
_symmetry_Int_Tables_number         221
_symmetry_cell_setting             cubic

loop_
_symmetry_equiv_pos_as_xyz
'+x,+y,+z'
'+z,+x,+y'
'+y,+z,+x'
'+x,+y,-z'
'+z,+x,-y'
'+y,+z,-x'
'-x,+y,+z'
'-z,+x,+y'
'-y,+z,+x'
'-x,+y,-z'
'-z,+x,-y'
'-y,+z,-x'
'+y,+x,+z'
'+x,+z,+y'
'+z,+y,+x'
'+y,+x,-z'
'+x,+z,-y'
'+z,+y,-x'
'+y,-x,+z'
'+x,-z,+y'
'+z,-y,+x'
'+y,-x,-z'
'+x,-z,-y'
'+z,-y,-x'
'-x,-y,-z'
'-z,-x,-y'
'-y,-z,-x'
'-x,-y,+z'
'-z,-x,+y'
'-y,-z,+x'
'+x,-y,-z'
'+z,-x,-y'
'+y,-z,-x'
'+x,-y,+z'
'+z,-x,+y'
'+y,-z,+x'
'-y,-x,-z'
'-x,-z,-y'
'-z,-y,-x'
'-y,-x,+z'
'-x,-z,+y'
'-z,-y,+x'
'-y,+x,-z'
'-x,+z,-y'
'-z,+y,-x'
'-y,+x,+z'
'-x,+z,+y'
'-z,+y,+x'

loop_
_atom_site_label
_atom_site_type_symbol
_atom_site_fract_x
_atom_site_fract_y
_atom_site_fract_z
    O1    O     0.0000    0.2122    0.5000
    O2    O     0.1103    0.1103    0.3384
    O3    O     0.0000    0.2967    0.2967
    T1    Si    0.0000    0.1823    0.3684
)CIF";

constexpr const char* kPtyCif = R"CIF(data_PTY

#**************************************************************************
#
# CIF taken from the IZA-SC Database of Zeolite Structures
# Ch. Baerlocher and L.B. McCusker
# Database of Zeolite Structures: http://www.iza-structure.org/databases/
#
# The atom coordinates and the cell parameters were optimized with DLS76
# assuming a pure SiO2 composition.
#
#**************************************************************************

_cell_length_a                   7.1158
_cell_length_b                   9.1748
_cell_length_c                   9.1160
_cell_angle_alpha               84.6450
_cell_angle_beta                83.8160
_cell_angle_gamma               86.6590

_symmetry_space_group_name_H-M     'P -1'
_symmetry_Int_Tables_number         2
_symmetry_cell_setting             triclinic

loop_
_symmetry_equiv_pos_as_xyz
'+x,+y,+z'
'-x,-y,-z'

loop_
_atom_site_label
_atom_site_type_symbol
_atom_site_fract_x
_atom_site_fract_y
_atom_site_fract_z
T1    Si    0.6596    0.3547    0.1505
T2    Si    0.2383    0.3911    0.0877
T3    Si    0.0289    0.8245    0.0166
T4    Si    0.0387    0.6195    0.2988
T5    Si    0.2851    0.5520    0.5500
O1    O     0.7179    0.3434 	0.3172
O2    O     0.7193    0.5098 	0.0663
O3    O     0.7674    0.2237 	0.0642
O4    O     0.1040    0.4831 	0.2032
O5    O     0.1345    0.2467 	0.0586
O6    O     0.0185    0.7665 	0.1894
O7    O     0.1626    0.4120 	0.6082
O8    O     0.1928    0.6392 	0.4115
O9    O     0.4343    0.3421 	0.1556
O10   O     0.0000    0.0000    0.0000
O11   O     0.5000    0.5000    0.5000
)CIF";
// clang-format on

} // namespace

TEST_CASE("readCif parses the real IZA LTA structure (cubic -- 48 symmetry operations)",
          "[validation][io][cif]") {
    const TempCifFile file("aleator_test_LTA.cif", kLtaCif);
    const StructureData structure = readCif(file.path());

    const auto& m = structure.lattice.matrix();
    REQUIRE(std::abs(vectorLength(m[0]) - 11.9190) < 1e-3);
    REQUIRE(std::abs(vectorLength(m[1]) - 11.9190) < 1e-3);
    REQUIRE(std::abs(vectorLength(m[2]) - 11.9190) < 1e-3);
    REQUIRE(std::abs(angleBetweenDegrees(m[0], m[1]) - 90.0) < 1e-3);
    REQUIRE(std::abs(angleBetweenDegrees(m[1], m[2]) - 90.0) < 1e-3);
    REQUIRE(std::abs(angleBetweenDegrees(m[0], m[2]) - 90.0) < 1e-3);

    const int siCount = countSpecies(structure, "Si");
    const int oCount = countSpecies(structure, "O");

    // Space-group-independent correctness proof: every zeolite framework is
    // a 4-connected net of TO4 tetrahedra, so O must be exactly 2x T.
    REQUIRE(oCount == 2 * siCount);

    // Published IZA framework composition for LTA in this (Pm-3m,
    // a=11.919 A) reduced cell: 24 T-atoms, 48 O per unit cell.
    REQUIRE(siCount == 24);
    REQUIRE(oCount == 48);
    REQUIRE(structure.particles.size() == 72);
}

TEST_CASE(
    "readCif parses the real IZA PTY structure (triclinic P-1) and expansion multiplicities "
    "match hand-derived site symmetry",
    "[validation][io][cif]") {
    // PTY's symmetry loop is just {identity, inversion through the origin}.
    // A generic asymmetric-unit atom (not itself at an inversion center)
    // has orbit size 2. Of the 16 listed sites, O10 = (0,0,0) and
    // O11 = (0.5,0.5,0.5) are each exactly their own inversion image
    // (-P == P, mod 1, since 2*0=0 and 2*0.5=1), so they sit ON inversion
    // centers and have orbit size 1; all 14 others (T1-T5, O1-O9) are
    // generic and have orbit size 2. Expected: T = 5*2 = 10,
    // O = 9*2 + 1 + 1 = 20, total = 30. (And 20 == 2*10, the same
    // topological check as the LTA case.) A parser that naively multiplied
    // every site by the symmetry-operation count without deduplicating
    // special positions would instead produce T=10, O=22 (double-counting
    // O10 and O11) — this test specifically catches that bug.
    const TempCifFile file("aleator_test_PTY.cif", kPtyCif);
    const StructureData structure = readCif(file.path());

    const auto& m = structure.lattice.matrix();
    REQUIRE(std::abs(vectorLength(m[0]) - 7.1158) < 1e-3);
    REQUIRE(std::abs(vectorLength(m[1]) - 9.1748) < 1e-3);
    REQUIRE(std::abs(vectorLength(m[2]) - 9.1160) < 1e-3);
    REQUIRE(std::abs(angleBetweenDegrees(m[1], m[2]) - 84.6450) < 1e-3);
    REQUIRE(std::abs(angleBetweenDegrees(m[0], m[2]) - 83.8160) < 1e-3);
    REQUIRE(std::abs(angleBetweenDegrees(m[0], m[1]) - 86.6590) < 1e-3);

    const int siCount = countSpecies(structure, "Si");
    const int oCount = countSpecies(structure, "O");

    REQUIRE(oCount == 2 * siCount);
    REQUIRE(siCount == 10);
    REQUIRE(oCount == 20);
    REQUIRE(structure.particles.size() == 30);

    // Occupancy and source-label bookkeeping should be complete and sane
    // even though this file has no _atom_site_occupancy column.
    REQUIRE(structure.occupancy.size() == structure.particles.size());
    for (double occ : structure.occupancy) {
        REQUIRE(occ == 1.0);
    }
    REQUIRE(structure.sourceLabel.size() == structure.particles.size());
}

TEST_CASE("readCif throws a clear error on a malformed CIF instead of crashing or mis-parsing",
          "[validation][io][cif]") {
    const std::string malformed = R"CIF(data_BROKEN
_cell_length_a    10.000
_cell_length_b    10.000
_cell_length_c    10.000
_cell_angle_alpha 90.0
_cell_angle_beta  90.0
_cell_angle_gamma 90.0

loop_
_atom_site_label
_atom_site_fract_x
_atom_site_fract_y
_atom_site_fract_z
Si1  0.0  0.0
)CIF";
    const TempCifFile file("aleator_test_broken.cif", malformed);
    REQUIRE_THROWS_AS(readCif(file.path()), aleator::io::CifParseError);
}

TEST_CASE("readCif throws a clear error when a required cell tag is missing",
          "[validation][io][cif]") {
    const std::string malformed = R"CIF(data_BROKEN
_cell_length_a    10.000
_cell_length_b    10.000
_cell_angle_alpha 90.0
_cell_angle_beta  90.0
_cell_angle_gamma 90.0

loop_
_atom_site_label
_atom_site_fract_x
_atom_site_fract_y
_atom_site_fract_z
Si1  0.0  0.0  0.0
)CIF";
    const TempCifFile file("aleator_test_missing_tag.cif", malformed);
    REQUIRE_THROWS_AS(readCif(file.path()), aleator::io::CifParseError);
}

TEST_CASE("readCif throws rather than opening a nonexistent file silently",
          "[validation][io][cif]") {
    const auto missing =
        std::filesystem::temp_directory_path() / "aleator_test_does_not_exist_12345.cif";
    REQUIRE_THROWS(readCif(missing));
}
