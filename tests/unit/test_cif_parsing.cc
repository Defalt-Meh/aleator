// io/ milestone: low-level CIF building blocks (tokenizer/document model,
// symmetry operation parsing, element table). End-to-end parsing of real
// structure files is validated in
// tests/validation/test_cif_real_structures.cc.

#include <catch2/catch_test_macros.hpp>

#include <cmath>
#include <stdexcept>
#include <string>

#include "io/cif_document.hpp"
#include "io/cif_error.hpp"
#include "io/element_table.hpp"
#include "io/symmetry_operation.hpp"

using aleator::io::CifDocument;
using aleator::io::CifParseError;
using aleator::io::parseCifDocument;
using aleator::io::parseCifNumber;
using aleator::io::tryParseCifNumber;

namespace {
constexpr double kEps = 1e-12;
}

TEST_CASE("parseCifNumber strips a trailing standard-uncertainty suffix", "[unit][io][cif]") {
    REQUIRE(std::abs(parseCifNumber("11.9190(0)") - 11.9190) < kEps);
    REQUIRE(std::abs(parseCifNumber("1.234(12)") - 1.234) < kEps);
    REQUIRE(std::abs(parseCifNumber("-0.5") - (-0.5)) < kEps);
}

TEST_CASE("tryParseCifNumber treats '?' and '.' as absent rather than erroring",
          "[unit][io][cif]") {
    REQUIRE_FALSE(tryParseCifNumber("?").has_value());
    REQUIRE_FALSE(tryParseCifNumber(".").has_value());
    REQUIRE(tryParseCifNumber("0.5").value() == 0.5);
}

TEST_CASE("parseCifNumber throws a clear error on garbage input", "[unit][io][cif]") {
    REQUIRE_THROWS_AS(parseCifNumber("abc"), CifParseError);
}

TEST_CASE("parseCifDocument reads single tags, quoted values, and a loop_", "[unit][io][cif]") {
    const std::string text = R"CIF(
data_TEST
_cell_length_a    10.000
_symmetry_space_group_name_H-M   'P 1'

loop_
_atom_site_label
_atom_site_fract_x
_atom_site_fract_y
_atom_site_fract_z
Si1  0.1  0.2  0.3
O1   0.4  0.5  0.6
)CIF";

    const CifDocument doc = parseCifDocument(text);
    REQUIRE(doc.blockName == "TEST");
    REQUIRE(doc.findTag("_cell_length_a").value() == "10.000");
    REQUIRE(doc.findTag("_symmetry_space_group_name_h-m").value() == "P 1");

    const auto* loop =
        doc.findLoopWithColumns({"_atom_site_label", "_atom_site_fract_x", "_atom_site_fract_y",
                                  "_atom_site_fract_z"});
    REQUIRE(loop != nullptr);
    REQUIRE(loop->rows.size() == 2);
    const auto labelCol = static_cast<std::size_t>(loop->columnIndex("_atom_site_label"));
    const auto zCol = static_cast<std::size_t>(loop->columnIndex("_atom_site_fract_z"));
    REQUIRE(loop->rows[0][labelCol] == "Si1");
    REQUIRE(loop->rows[1][zCol] == "0.6");
}

TEST_CASE("parseCifDocument handles tabs as column separators", "[unit][io][cif]") {
    const std::string text =
        "data_TEST\n_cell_length_a 1.0\nloop_\n_atom_site_label\n_atom_site_fract_x\nO1\t0.5\n";
    const CifDocument doc = parseCifDocument(text);
    const auto* loop = doc.findLoopWithColumns({"_atom_site_label", "_atom_site_fract_x"});
    REQUIRE(loop != nullptr);
    REQUIRE(loop->rows.at(0).at(1) == "0.5");
}

TEST_CASE("parseCifDocument throws on a loop_ row that doesn't divide evenly into its columns",
          "[unit][io][cif]") {
    const std::string text = R"CIF(
data_TEST
loop_
_atom_site_label
_atom_site_fract_x
_atom_site_fract_y
_atom_site_fract_z
Si1  0.1  0.2
)CIF";
    REQUIRE_THROWS_AS(parseCifDocument(text), CifParseError);
}

TEST_CASE("parseCifDocument throws on an unterminated quoted string", "[unit][io][cif]") {
    const std::string text = "data_TEST\n_tag   'unterminated\n";
    REQUIRE_THROWS_AS(parseCifDocument(text), CifParseError);
}

TEST_CASE("parseCifDocument throws when no data_ block is present", "[unit][io][cif]") {
    const std::string text = "_cell_length_a 10.0\n";
    REQUIRE_THROWS_AS(parseCifDocument(text), CifParseError);
}

TEST_CASE("parseCifDocument throws on unrecognized top-level content", "[unit][io][cif]") {
    const std::string text = "data_TEST\nthis is not valid CIF syntax at all\n";
    REQUIRE_THROWS_AS(parseCifDocument(text), CifParseError);
}

TEST_CASE("parseSymmetryOperation: identity", "[unit][io][cif]") {
    const auto op = aleator::io::parseSymmetryOperation("x,y,z");
    const auto result = op.apply({0.2, 0.3, 0.4});
    REQUIRE(std::abs(result[0] - 0.2) < kEps);
    REQUIRE(std::abs(result[1] - 0.3) < kEps);
    REQUIRE(std::abs(result[2] - 0.4) < kEps);
}

TEST_CASE("parseSymmetryOperation: inversion", "[unit][io][cif]") {
    const auto op = aleator::io::parseSymmetryOperation("-x,-y,-z");
    const auto result = op.apply({0.2, 0.3, 0.4});
    REQUIRE(std::abs(result[0] + 0.2) < kEps);
    REQUIRE(std::abs(result[1] + 0.3) < kEps);
    REQUIRE(std::abs(result[2] + 0.4) < kEps);
}

TEST_CASE("parseSymmetryOperation: fractional translations", "[unit][io][cif]") {
    const auto op = aleator::io::parseSymmetryOperation("1/2+x,1/2-y,z");
    const auto result = op.apply({0.1, 0.1, 0.1});
    REQUIRE(std::abs(result[0] - 0.6) < kEps);  // 0.5 + 0.1
    REQUIRE(std::abs(result[1] - 0.4) < kEps);  // 0.5 - 0.1
    REQUIRE(std::abs(result[2] - 0.1) < kEps);
}

TEST_CASE("parseSymmetryOperation: trigonal-style combined coefficients", "[unit][io][cif]") {
    // "-y,x-y,z+1/3" is a real hexagonal/trigonal-group symmetry operator.
    const auto op = aleator::io::parseSymmetryOperation("-y,x-y,z+1/3");
    const auto result = op.apply({0.5, 0.2, 0.7});
    REQUIRE(std::abs(result[0] - (-0.2)) < kEps);        // -y
    REQUIRE(std::abs(result[1] - (0.5 - 0.2)) < kEps);   // x - y
    REQUIRE(std::abs(result[2] - (0.7 + 1.0 / 3.0)) < kEps); // z + 1/3
}

TEST_CASE("parseSymmetryOperation throws on malformed input", "[unit][io][cif]") {
    REQUIRE_THROWS_AS(aleator::io::parseSymmetryOperation("x,y"), CifParseError);
    REQUIRE_THROWS_AS(aleator::io::parseSymmetryOperation("x,y,w"), CifParseError);
}

TEST_CASE("atomicMass returns known values and nullopt for unrecognized symbols",
          "[unit][io][cif]") {
    REQUIRE(aleator::io::atomicMass("Si").has_value());
    REQUIRE(std::abs(*aleator::io::atomicMass("Si") - 28.085) < 1e-6);
    REQUIRE(aleator::io::atomicMass("O").has_value());
    REQUIRE_FALSE(aleator::io::atomicMass("Xx").has_value());
}

TEST_CASE("elementFromLabel strips numeric/charge suffixes and disambiguates 1- vs 2-letter symbols",
          "[unit][io][cif]") {
    REQUIRE(aleator::io::elementFromLabel("Si1").value() == "Si");
    REQUIRE(aleator::io::elementFromLabel("O12A").value() == "O");
    REQUIRE(aleator::io::elementFromLabel("Cu2+").value() == "Cu");
    REQUIRE(aleator::io::elementFromLabel("Cl1-").value() == "Cl");
    REQUIRE_FALSE(aleator::io::elementFromLabel("1Si").has_value());
}
