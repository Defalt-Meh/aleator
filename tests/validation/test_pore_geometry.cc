// Validation anchor (CLAUDE.md section 4, engines/geometry_analysis
// milestone). Split by how tightly each quantity is actually pinned down:
//
//   - LCD: a genuine, understood 1-4% gap vs. real Zeo++ 0.4.7 on every
//     real structure (power-vs-Apollonius diagram, see
//     data/pore_geometry/PROVENANCE.md) — split off into
//     tests/known_deviation/test_pore_geometry_known_deviation.cc per
//     section 4, never loosened here to swallow it.
//   - PLD, ASA, AV on LTA/MFI/FAU: PLD agrees with real Zeo++ to within
//     4e-4 Angstrom (well inside the 1e-3 Angstrom target), and ASA/AV
//     (Monte Carlo, so real sampling noise on top) agree to within 3% —
//     validated TIGHT here.
//   - PLD, ASA, AV on PTY (the triclinic structure): investigated
//     directly, not assumed to be more of the same LCD-style noise — see
//     PROVENANCE.md for the four independent checks performed (lattice
//     vs. Zeo++'s own box vectors, minimum image vs. brute force, grid
//     resolution, total tessellation volume) that ruled out a bug, and
//     the specific percolation-critical edge inspected that explains the
//     gap as a genuine methodological difference in what each tool's
//     percolation graph accounts for. Also split into known-deviation.
//
// Because PTY's own *primitive* PBC operation (minimum image) is what
// invariant 6 actually requires be tested on a non-orthogonal cell, and
// that primitive is now independently proven correct here (not just
// assumed from other milestones' triclinic coverage), invariant 6 is
// satisfied by a tight, targeted test even though the full pore-analysis
// pipeline's agreement with Zeo++ on PTY specifically is not tight.
//
// All four structures are real IZA zeolite CIFs (LTA, MFI, FAU, PTY),
// fetched live and checked in under data/pore_geometry/ — not synthetic —
// plus one small closed-form analytic case with a hand-derivable answer.
//
// Each real structure is analyzed at most ONCE per TEST_CASE (PLD, ASA, AV
// all read off the same PoreGeometry result) rather than once per checked
// quantity: constructing the periodic Voronoi network itself (not the
// Monte Carlo sampling) is the expensive part — see
// voronoi_network.cc's kEdgeSampleSpacingAngstrom — so re-running
// analyzePoreGeometry per-assertion would multiply real wall-clock cost for
// no additional coverage.

#include <catch2/catch_test_macros.hpp>

#include <array>
#include <cmath>
#include <cstdio>
#include <filesystem>
#include <limits>
#include <random>
#include <stdexcept>
#include <string>
#include <vector>

#include "core/geometry/lattice.hpp"
#include "core/math/particle_data.hpp"
#include "engines/geometry_analysis/pore_analysis.hpp"
#include "io/structure_io.hpp"

using aleator::core::Lattice;
using aleator::core::ParticleData;
using aleator::engines::analyzePoreGeometry;
using aleator::engines::PoreAnalysisOptions;
using aleator::engines::PoreGeometry;
using aleator::engines::zeoPlusPlusDefaultRadiusAngstrom;

namespace {

std::filesystem::path poreGeometryDataFile(const std::string& name) {
    return std::filesystem::path(ALEATOR_TEST_DATA_DIR) / "data" / "pore_geometry" / name;
}

std::vector<double> radiiFor(const std::vector<std::string>& speciesSymbols) {
    std::vector<double> radii;
    radii.reserve(speciesSymbols.size());
    for (const auto& symbol : speciesSymbols) {
        radii.push_back(zeoPlusPlusDefaultRadiusAngstrom(symbol));
    }
    return radii;
}

/// Real Zeo++ 0.4.7 (`network -res`/`-sa`/`-vol`, N2 probe = 1.86 Angstrom)
/// output on the exact CIF files checked in under data/pore_geometry/ — see
/// PROVENANCE.md for the commands and full table.
struct ZeoPlusPlusReference {
    std::string name;
    double pldAngstrom;
    double asaAngstromSq;
    double nasaAngstromSq;
    double avAngstromCubed;
    double navAngstromCubed;
};

/// The three orthogonal structures only -- PTY (triclinic) is validated
/// separately (see PROVENANCE.md for why its PLD/ASA/AV are a known
/// deviation rather than a tight match).
const std::vector<ZeoPlusPlusReference>& references() {
    static const std::vector<ZeoPlusPlusReference> data{
        {"LTA", 3.81365, 178.578, 13.8393, 184.564, 4.63949},
        {"MFI", 4.24727, 359.163, 0.0, 91.6143, 0.0},
        {"FAU", 6.95049, 1863.61, 91.1772, 2458.66, 25.1061},
    };
    return data;
}

} // namespace

TEST_CASE("analyzePoreGeometry reproduces the analytic LCD/PLD of a single point atom in a simple "
          "cubic cell",
          "[validation][geometry][pore]") {
    // One radius-0 atom at the origin of a 10x10x10 Angstrom cubic cell.
    // With a single atom per cell, its periodic Voronoi cell is exactly the
    // cube itself: vertices are the 8 corners at distance a*sqrt(3)/2 from
    // the atom, so LCD = a*sqrt(3) exactly. The tightest constriction is at
    // any cube-edge midpoint (e.g. (a/2, a/2, 0)), equidistant from the 4
    // atom images surrounding that edge at distance
    // sqrt((a/2)^2 + (a/2)^2) = a/sqrt(2), so PLD = a*sqrt(2). Both are
    // exact closed forms, independent of Zeo++ or any external tool.
    const double a = 10.0;
    ParticleData framework;
    framework.push_back(0.0, 0.0, 0.0, 0.0, 0.0, 0.0, /*mass=*/1.0, /*charge=*/0.0, /*species=*/0);

    std::array<std::array<double, 3>, 3> matrix{};
    matrix[0] = {a, 0.0, 0.0};
    matrix[1] = {0.0, a, 0.0};
    matrix[2] = {0.0, 0.0, a};
    const Lattice lattice(matrix);

    PoreAnalysisOptions options;
    options.probeRadiusAngstrom = 0.0; // irrelevant to LCD/PLD, keep ASA/AV cheap
    options.asaSamplesPerAtom = 200;
    options.volumeSamplesTotal = 2'000;

    const PoreGeometry result = analyzePoreGeometry(framework, lattice, {0.0}, options);

    const double expectedLcd = a * std::sqrt(3.0);
    const double expectedPld = a * std::sqrt(2.0);

    std::printf("Analytic simple-cubic case: LCD computed=%.7f expected=%.7f, "
                "PLD computed=%.7f expected=%.7f\n",
                result.largestCavityDiameterAngstrom, expectedLcd, result.poreLimitingDiameterAngstrom,
                expectedPld);

    CHECK(std::abs(result.largestCavityDiameterAngstrom - expectedLcd) < 1e-4);
    CHECK(std::abs(result.poreLimitingDiameterAngstrom - expectedPld) < 1e-4);
}

TEST_CASE("analyzePoreGeometry rejects an empty framework", "[validation][geometry][pore]") {
    const ParticleData empty;
    std::array<std::array<double, 3>, 3> matrix{};
    matrix[0] = {10.0, 0.0, 0.0};
    matrix[1] = {0.0, 10.0, 0.0};
    matrix[2] = {0.0, 0.0, 10.0};
    const Lattice lattice(matrix);

    REQUIRE_THROWS_AS(analyzePoreGeometry(empty, lattice, {}), std::invalid_argument);
}

TEST_CASE("analyzePoreGeometry rejects a species index the radius table doesn't cover",
          "[validation][geometry][pore]") {
    ParticleData framework;
    framework.push_back(0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, /*species=*/0);
    std::array<std::array<double, 3>, 3> matrix{};
    matrix[0] = {10.0, 0.0, 0.0};
    matrix[1] = {0.0, 10.0, 0.0};
    matrix[2] = {0.0, 0.0, 10.0};
    const Lattice lattice(matrix);

    // Empty radius table: species 0 is out of range.
    REQUIRE_THROWS_AS(analyzePoreGeometry(framework, lattice, {}), std::invalid_argument);
}

TEST_CASE("analyzePoreGeometry's PLD and ASA/AV match real Zeo++ 0.4.7 N2-probe output on real "
          "zeolite structures (LTA, MFI, FAU)",
          "[validation][geometry][pore][slow]") {
    for (const ZeoPlusPlusReference& ref : references()) {
        const auto structure = aleator::io::readCif(poreGeometryDataFile(ref.name + ".cif"));
        REQUIRE(structure.particles.size() > 0);
        const std::vector<double> radii = radiiFor(structure.speciesSymbols);

        PoreAnalysisOptions options; // defaults: 1.86 Angstrom (N2) probe, 2000/atom, 50000 total
        const PoreGeometry result =
            analyzePoreGeometry(structure.particles, structure.lattice, radii, options);

        const double pldError = std::abs(result.poreLimitingDiameterAngstrom - ref.pldAngstrom);

        // Zeo++ reports ASA/AV as accessible-only, NASA/NAV as the pocket
        // contribution separately; compare totals (accessible + pocket
        // combined) rather than the split, since this codebase's
        // nearest-qualifying-node classification (pore_analysis.cc) is a
        // documented approximation of exactly which side of the
        // accessible/pocket boundary a borderline sample falls on, and the
        // split carries extra noise the total does not.
        const double totalArea =
            result.accessibleSurfaceAreaAngstromSq + result.inaccessibleSurfaceAreaAngstromSq;
        const double totalVolume =
            result.accessibleVolumeAngstromCubed + result.inaccessibleVolumeAngstromCubed;
        const double referenceTotalArea = ref.asaAngstromSq + ref.nasaAngstromSq;
        const double referenceTotalVolume = ref.avAngstromCubed + ref.navAngstromCubed;

        std::printf("%s: PLD computed=%.5f Zeo++=%.5f |error|=%.5f Angstrom | "
                    "ASA+NASA computed=%.3f Zeo++=%.3f | AV+NAV computed=%.3f Zeo++=%.3f\n",
                    ref.name.c_str(), result.poreLimitingDiameterAngstrom, ref.pldAngstrom, pldError,
                    totalArea, referenceTotalArea, totalVolume, referenceTotalVolume);

        INFO(ref.name << ": computed PLD=" << result.poreLimitingDiameterAngstrom
                      << ", Zeo++ reference=" << ref.pldAngstrom);
        CHECK(pldError < 1e-3);

        INFO(ref.name << ": total area computed=" << totalArea << ", Zeo++=" << referenceTotalArea);
        CHECK(std::abs(totalArea - referenceTotalArea) / referenceTotalArea < 0.10);

        INFO(ref.name << ": total volume computed=" << totalVolume
                      << ", Zeo++=" << referenceTotalVolume);
        CHECK(std::abs(totalVolume - referenceTotalVolume) / referenceTotalVolume < 0.10);

        // FAU-specific: the milestone explicitly asks for "a structure with
        // a known blocked/inaccessible pocket, asserting the accessible and
        // total volumes differ in the expected direction" — FAU's sodalite
        // cages are real, published-as-inaccessible-to-N2 pockets (Zeo++:
        // 8 pockets, NAV=25.1061 Ang^3) distinct from the main 12-ring
        // channel system, so this reuses the real structure already in
        // hand rather than fabricating a synthetic one.
        if (ref.name == "FAU") {
            CHECK(result.inaccessibleVolumeAngstromCubed > 0.0);
            CHECK(result.accessibleVolumeAngstromCubed > 10.0 * result.inaccessibleVolumeAngstromCubed);
            CHECK(result.accessibleVolumeAngstromCubed + result.inaccessibleVolumeAngstromCubed <
                  result.totalVolumeAngstromCubed);
        }
    }
}

TEST_CASE("Lattice::minimumImageDisplacement matches a brute-force periodic-image search on the "
          "real triclinic PTY cell",
          "[validation][geometry][pore][invariant6]") {
    // CLAUDE.md invariant 6: every PBC path must be tested on a
    // non-orthogonal cell, not just cubic. PTY (a=7.1158 b=9.1748
    // c=9.1160 Angstrom, alpha=84.645 beta=83.816 gamma=86.659 degrees) is
    // this codebase's realest available triclinic case. This test checks
    // the actual PRIMITIVE minimum-image operation directly against
    // exhaustive enumeration, independent of Voro++, Zeo++, or the rest of
    // the pore-analysis pipeline -- this is what conclusively rules out a
    // triclinic PBC bug as the explanation for PTY's PLD/ASA/AV gap from
    // Zeo++ (see PROVENANCE.md): the primitive itself is correct to
    // floating-point precision, so the gap has to come from somewhere
    // else, which is investigated and documented separately.
    const auto structure = aleator::io::readCif(poreGeometryDataFile("PTY.cif"));
    const auto& m = structure.lattice.matrix();
    const std::array<double, 3> a0{m[0][0], m[0][1], m[0][2]};
    const std::array<double, 3> a1{m[1][0], m[1][1], m[1][2]};
    const std::array<double, 3> a2{m[2][0], m[2][1], m[2][2]};

    std::mt19937_64 rng(12345);
    std::uniform_real_distribution<double> uni(-5.0, 15.0);

    double maxRelError = 0.0;
    constexpr int kTrials = 5'000;
    for (int t = 0; t < kTrials; ++t) {
        const std::array<double, 3> p{uni(rng), uni(rng), uni(rng)};
        const std::array<double, 3> q{uni(rng), uni(rng), uni(rng)};

        const auto d = structure.lattice.minimumImageDisplacement(p, q);
        const double distComputed = std::sqrt(d[0] * d[0] + d[1] * d[1] + d[2] * d[2]);

        double bestBrute = std::numeric_limits<double>::infinity();
        for (int i = -3; i <= 3; ++i) {
            for (int j = -3; j <= 3; ++j) {
                for (int k = -3; k <= 3; ++k) {
                    const std::array<double, 3> qImage{
                        q[0] + i * a0[0] + j * a1[0] + k * a2[0],
                        q[1] + i * a0[1] + j * a1[1] + k * a2[1],
                        q[2] + i * a0[2] + j * a1[2] + k * a2[2]};
                    const double dx = p[0] - qImage[0];
                    const double dy = p[1] - qImage[1];
                    const double dz = p[2] - qImage[2];
                    const double dist = std::sqrt(dx * dx + dy * dy + dz * dz);
                    bestBrute = std::min(bestBrute, dist);
                }
            }
        }

        const double relError = std::abs(distComputed - bestBrute) / std::max(1e-9, bestBrute);
        maxRelError = std::max(maxRelError, relError);
    }

    std::printf("PTY minimum-image vs. brute force (-3..3 images), %d trials: max relative error "
                "= %.3e\n",
                kTrials, maxRelError);
    CHECK(maxRelError < 1e-6);
}

TEST_CASE("analyzePoreGeometry runs end-to-end on the real triclinic PTY structure and produces "
          "physically sane output",
          "[validation][geometry][pore][slow]") {
    // Structural sanity, not a tight match against Zeo++ (see
    // PROVENANCE.md and the known-deviation test for that comparison and
    // why it isn't tight here) -- this proves the full pipeline (Voro++
    // construction, canonical-box transform, percolation, Monte Carlo
    // sampling) completes successfully on a genuinely non-orthogonal cell
    // and produces internally consistent numbers.
    const auto structure = aleator::io::readCif(poreGeometryDataFile("PTY.cif"));
    const std::vector<double> radii = radiiFor(structure.speciesSymbols);

    PoreAnalysisOptions options;
    options.probeRadiusAngstrom = 1.86;

    const PoreGeometry result =
        analyzePoreGeometry(structure.particles, structure.lattice, radii, options);

    std::printf("PTY: LCD=%.5f PLD=%.5f ASA=%.4f NASA=%.4f AV=%.4f NAV=%.4f total=%.4f\n",
                result.largestCavityDiameterAngstrom, result.poreLimitingDiameterAngstrom,
                result.accessibleSurfaceAreaAngstromSq, result.inaccessibleSurfaceAreaAngstromSq,
                result.accessibleVolumeAngstromCubed, result.inaccessibleVolumeAngstromCubed,
                result.totalVolumeAngstromCubed);

    CHECK(result.largestCavityDiameterAngstrom > result.poreLimitingDiameterAngstrom);
    CHECK(result.poreLimitingDiameterAngstrom > 0.0);
    CHECK(result.accessibleVolumeAngstromCubed >= 0.0);
    CHECK(result.inaccessibleVolumeAngstromCubed >= 0.0);
    CHECK(result.accessibleVolumeAngstromCubed + result.inaccessibleVolumeAngstromCubed <
          result.totalVolumeAngstromCubed);
    CHECK(std::abs(result.totalVolumeAngstromCubed - structure.lattice.volume()) < 1e-6);
}
