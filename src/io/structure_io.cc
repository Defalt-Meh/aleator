#include "io/structure_io.hpp"

#include <array>
#include <cmath>
#include <cstdint>
#include <fstream>
#include <numbers>
#include <sstream>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

#include "core/exceptions.hpp"
#include "io/cif_document.hpp"
#include "io/cif_error.hpp"
#include "io/element_table.hpp"
#include "io/symmetry_operation.hpp"

namespace aleator::io {

namespace {

using Vec3 = std::array<double, 3>;

/// Fractional-to-fractional dedup tolerance for symmetry-orbit expansion,
/// in Å (Cartesian, after converting the wrapped fractional difference
/// through the cell matrix). Real CIF coordinates are typically given to
/// 4-5 decimal fractional digits, so 1e-3 Å comfortably separates distinct
/// atoms from rounding noise on a shared special position.
constexpr double kOrbitDedupToleranceAngstrom = 1e-3;

std::string readWholeFile(const std::filesystem::path& file) {
    std::ifstream in(file, std::ios::binary);
    if (!in) {
        throw CifParseError("could not open file: " + file.string());
    }
    std::ostringstream buffer;
    buffer << in.rdbuf();
    return buffer.str();
}

double requireCellNumber(const CifDocument& doc, const std::string& tag) {
    const auto value = doc.findTag(tag);
    if (!value.has_value()) {
        throw CifParseError("missing required tag " + tag);
    }
    return parseCifNumber(*value);
}

/// Builds the lattice matrix from the standard crystallographic a, b, c,
/// alpha, beta, gamma parameterization (CLAUDE.md #2.4: CIF angles are in
/// degrees, converted to radians here at the I/O boundary; lengths are
/// already in Å, matching Aleator's internal length unit, so no length
/// conversion is needed).
core::Lattice buildLattice(const CifDocument& doc) {
    const double a = requireCellNumber(doc, "_cell_length_a");
    const double b = requireCellNumber(doc, "_cell_length_b");
    const double c = requireCellNumber(doc, "_cell_length_c");
    const double alphaDeg = requireCellNumber(doc, "_cell_angle_alpha");
    const double betaDeg = requireCellNumber(doc, "_cell_angle_beta");
    const double gammaDeg = requireCellNumber(doc, "_cell_angle_gamma");

    for (const auto& [name, value] :
         {std::pair{"a", a}, std::pair{"b", b}, std::pair{"c", c}}) {
        if (!(value > 0.0)) {
            throw CifParseError(std::string("_cell_length_") + name + " must be positive, got " +
                                 std::to_string(value));
        }
    }
    for (const auto& [name, value] : {std::pair{"alpha", alphaDeg}, std::pair{"beta", betaDeg},
                                       std::pair{"gamma", gammaDeg}}) {
        if (!(value > 0.0 && value < 180.0)) {
            throw CifParseError(std::string("_cell_angle_") + name +
                                 " must be strictly between 0 and 180 degrees, got " +
                                 std::to_string(value));
        }
    }

    constexpr double kDegToRad = std::numbers::pi / 180.0;
    const double alpha = alphaDeg * kDegToRad;
    const double beta = betaDeg * kDegToRad;
    const double gamma = gammaDeg * kDegToRad;
    const double cosA = std::cos(alpha);
    const double cosB = std::cos(beta);
    const double cosG = std::cos(gamma);
    const double sinG = std::sin(gamma);

    const Vec3 aVec{a, 0.0, 0.0};
    const Vec3 bVec{b * cosG, b * sinG, 0.0};

    const double cx = c * cosB;
    const double cy = c * (cosA - cosB * cosG) / sinG;
    const double czSquared = 1.0 - cosA * cosA - cosB * cosB - cosG * cosG + 2.0 * cosA * cosB * cosG;
    if (czSquared < 0.0) {
        throw CifParseError(
            "cell angles alpha/beta/gamma do not describe a valid parallelepiped "
            "(negative squared c-component)");
    }
    const Vec3 cVec{cx, cy, c * std::sqrt(czSquared) / sinG};

    return core::Lattice({{aVec, bVec, cVec}});
}

std::vector<SymmetryOperation> readSymmetryOperations(const CifDocument& doc) {
    const CifLoop* loop = doc.findLoopWithColumns({"_symmetry_equiv_pos_as_xyz"});
    std::string columnTag = "_symmetry_equiv_pos_as_xyz";
    if (loop == nullptr) {
        loop = doc.findLoopWithColumns({"_space_group_symop_operation_xyz"});
        columnTag = "_space_group_symop_operation_xyz";
    }
    if (loop == nullptr) {
        // No explicit symmetry loop: treat as P1 (identity only), which is
        // both the CIF default and how many already-fully-expanded MOF
        // CIFs are written.
        return {parseSymmetryOperation("x,y,z")};
    }

    const int column = loop->columnIndex(columnTag);
    std::vector<SymmetryOperation> ops;
    ops.reserve(loop->rows.size());
    for (const auto& row : loop->rows) {
        ops.push_back(parseSymmetryOperation(row[static_cast<std::size_t>(column)]));
    }
    if (ops.empty()) {
        throw CifParseError("symmetry loop (" + columnTag + ") contains no operations");
    }
    return ops;
}

Vec3 wrapToUnitCell(Vec3 f) {
    for (double& x : f) {
        x -= std::floor(x);
    }
    return f;
}

/// True if `candidate` coincides (within kOrbitDedupToleranceAngstrom) with
/// any position already in `existing`, accounting for periodic wraparound
/// (e.g. fractional 0.0 and 0.99999... are the same point).
bool coincidesWithAny(const core::Lattice& lattice, const Vec3& candidate,
                       const std::vector<Vec3>& existing) {
    for (const auto& other : existing) {
        Vec3 diff{candidate[0] - other[0], candidate[1] - other[1], candidate[2] - other[2]};
        for (double& d : diff) {
            d -= std::round(d);
        }
        const auto cartesianDiff = lattice.fractionalToCartesian(diff);
        const double distance = std::sqrt(cartesianDiff[0] * cartesianDiff[0] +
                                           cartesianDiff[1] * cartesianDiff[1] +
                                           cartesianDiff[2] * cartesianDiff[2]);
        if (distance < kOrbitDedupToleranceAngstrom) {
            return true;
        }
    }
    return false;
}

std::string requireElement(const std::string& label, const std::string* typeSymbol) {
    if (typeSymbol != nullptr) {
        if (const auto element = elementFromLabel(*typeSymbol)) {
            return *element;
        }
    }
    if (const auto element = elementFromLabel(label)) {
        return *element;
    }
    throw CifParseError("could not determine a chemical element for atom site \"" + label + "\"");
}

} // namespace

StructureData readCif(const std::filesystem::path& file) {
    const std::string text = readWholeFile(file);
    const CifDocument doc = parseCifDocument(text);

    StructureData result;
    result.lattice = buildLattice(doc);
    const std::vector<SymmetryOperation> symmetryOps = readSymmetryOperations(doc);

    const CifLoop* atomLoop = doc.findLoopWithColumns(
        {"_atom_site_label", "_atom_site_fract_x", "_atom_site_fract_y", "_atom_site_fract_z"});
    if (atomLoop == nullptr) {
        throw CifParseError(
            "no _atom_site loop with a label and fractional x/y/z coordinates was found");
    }

    const int labelCol = atomLoop->columnIndex("_atom_site_label");
    const int typeCol = atomLoop->columnIndex("_atom_site_type_symbol");
    const int xCol = atomLoop->columnIndex("_atom_site_fract_x");
    const int yCol = atomLoop->columnIndex("_atom_site_fract_y");
    const int zCol = atomLoop->columnIndex("_atom_site_fract_z");
    const int occCol = atomLoop->columnIndex("_atom_site_occupancy");
    const int chargeCol = atomLoop->columnIndex("_atom_site_charge");

    std::unordered_map<std::string, std::uint32_t> speciesIndex;

    for (const auto& row : atomLoop->rows) {
        const std::string& label = row[static_cast<std::size_t>(labelCol)];
        const std::string* typeSymbol =
            typeCol >= 0 ? &row[static_cast<std::size_t>(typeCol)] : nullptr;
        const std::string element = requireElement(label, typeSymbol);

        const auto mass = atomicMass(element);
        if (!mass.has_value()) {
            throw CifParseError("unknown element \"" + element +
                                 "\" (no atomic mass data) for atom site \"" + label + "\"");
        }

        const double fx = parseCifNumber(row[static_cast<std::size_t>(xCol)]);
        const double fy = parseCifNumber(row[static_cast<std::size_t>(yCol)]);
        const double fz = parseCifNumber(row[static_cast<std::size_t>(zCol)]);
        const double occupancy =
            occCol >= 0 ? tryParseCifNumber(row[static_cast<std::size_t>(occCol)]).value_or(1.0)
                        : 1.0;
        // Optional -- absent for the (already-validated) uncharged
        // methane/IRMOF-1 case, present for DDEC/DFT-partitioned partial
        // charges (e.g. the CO2/IRMOF-1 CRAFTED reference). Defaults to
        // 0.0, matching this reader's previous unconditional behavior.
        const double charge =
            chargeCol >= 0 ? tryParseCifNumber(row[static_cast<std::size_t>(chargeCol)]).value_or(0.0)
                            : 0.0;

        const auto speciesIt = speciesIndex.find(element);
        std::uint32_t speciesId{};
        if (speciesIt != speciesIndex.end()) {
            speciesId = speciesIt->second;
        } else {
            speciesId = static_cast<std::uint32_t>(result.speciesSymbols.size());
            speciesIndex.emplace(element, speciesId);
            result.speciesSymbols.push_back(element);
        }

        std::vector<Vec3> orbit;
        for (const auto& op : symmetryOps) {
            const Vec3 image = wrapToUnitCell(op.apply({fx, fy, fz}));
            if (!coincidesWithAny(result.lattice, image, orbit)) {
                orbit.push_back(image);
            }
        }

        for (const auto& fractional : orbit) {
            const auto cartesian = result.lattice.fractionalToCartesian(
                {fractional[0], fractional[1], fractional[2]});
            result.particles.push_back(cartesian[0], cartesian[1], cartesian[2], 0.0, 0.0, 0.0,
                                        *mass, charge, speciesId);
            result.occupancy.push_back(occupancy);
            result.sourceLabel.push_back(label);
        }
    }

    return result;
}

void writePdb(const std::filesystem::path&, const core::ParticleData&, const core::Lattice&) {
    throw aleator::NotImplemented("writePdb");
}

void writeLammpsData(const std::filesystem::path&, const core::ParticleData&,
                      const core::Lattice&) {
    throw aleator::NotImplemented("writeLammpsData");
}

void appendTrajectoryFrame(const std::filesystem::path&, const core::ParticleData&,
                            const core::Lattice&) {
    throw aleator::NotImplemented("appendTrajectoryFrame");
}

} // namespace aleator::io
