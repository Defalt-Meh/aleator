// Informational known-deviation test (CLAUDE.md section 4) for
// engines/geometry_analysis. Two distinct, separately-investigated gaps
// against real Zeo++ 0.4.7 output (see
// tests/validation/data/pore_geometry/PROVENANCE.md for the full
// writeup of both):
//
//   - LCD (largest cavity diameter), all four structures: a systematic
//     1-4% deviation from Voro++'s power (Laguerre) diagram vs. the
//     "true" Apollonius (additively-weighted) diagram most people mean by
//     "Voronoi decomposition of unequal spheres" -- not implemented by
//     any fast library, including the one Zeo++ itself uses. Ruled out as
//     a missed-candidate-point bug (edge-interior sampling, already
//     computed for PLD, never exceeded the vertex-only maximum on LTA).
//   - PLD/ASA/AV, PTY only: a *different*, specifically investigated
//     cause -- a real percolation-critical Voronoi edge whose interior
//     (sampled by this codebase's edge machinery, validated on LTA) dips
//     below both its endpoint vertices' radii, while Zeo++'s reported
//     threshold exactly matches one of this codebase's own vertex radii,
//     suggesting a methodological difference in what each tool's
//     percolation graph accounts for -- not a bug. Four independent
//     checks (lattice vs. Zeo++'s own box vectors, minimum image vs.
//     brute force over 200,000 trials, Voro++ grid-resolution
//     insensitivity, exact total-tessellation-volume match) ruled out a
//     bug in this codebase's own geometry/PBC code before this was
//     accepted as a known deviation rather than debugged further.
//
// PLD/ASA/AV on LTA/MFI/FAU are NOT here: they're validated TIGHT in
// tests/validation/test_pore_geometry.cc (agree with Zeo++ to <=4e-4
// Angstrom / <3%). This file exists so a genuine future regression in
// either gap is still caught, without the tight suite being loosened to
// swallow either one.

#include <catch2/catch_test_macros.hpp>

#include <cmath>
#include <cstdio>
#include <filesystem>
#include <string>
#include <vector>

#include "engines/geometry_analysis/pore_analysis.hpp"
#include "io/structure_io.hpp"

using aleator::engines::analyzePoreGeometry;
using aleator::engines::PoreAnalysisOptions;
using aleator::engines::PoreGeometry;
using aleator::engines::zeoPlusPlusDefaultRadiusAngstrom;

namespace {

std::filesystem::path poreGeometryDataFile(const std::string& name) {
    return std::filesystem::path(ALEATOR_KNOWN_DEVIATION_DATA_DIR) / "data" / "pore_geometry" / name;
}

std::vector<double> radiiFor(const std::vector<std::string>& speciesSymbols) {
    std::vector<double> radii;
    radii.reserve(speciesSymbols.size());
    for (const auto& symbol : speciesSymbols) {
        radii.push_back(zeoPlusPlusDefaultRadiusAngstrom(symbol));
    }
    return radii;
}

struct LcdBaseline {
    std::string name;
    double zeoPlusPlusLcdAngstrom;
    /// Recorded from a real run of this exact test (see
    /// known_deviation_baseline.md for provenance). Deterministic (no
    /// Monte Carlo involved in LCD), so no noise margin is needed on top.
    double baselineRelativeDeviation;
};

const std::vector<LcdBaseline>& lcdBaselines() {
    static const std::vector<LcdBaseline> data{
        {"LTA", 10.65071, 0.0362},
        {"MFI", 5.95565, 0.0044},
        {"FAU", 10.83869, 0.0126},
        {"PTY", 4.71404, 0.0080},
    };
    return data;
}

// Small additive margin (percentage points) on top of each baseline, since
// LCD/PLD are exact graph computations with no Monte Carlo noise but are
// not literally bitwise pinned to a specific Voro++ internal ordering
// across platforms/compilers -- CLAUDE.md section 4: sized to catch a
// genuine regression, not to swallow one.
constexpr double kAllowedDeviationGrowth = 0.02;

} // namespace

TEST_CASE("Pore geometry LCD's deviation from real Zeo++ 0.4.7 has not grown beyond its recorded "
          "baseline (LTA, MFI, FAU, PTY)",
          "[known-deviation][geometry][pore][slow]") {
    std::printf("%10s %14s %14s %14s %14s\n", "structure", "computed", "Zeo++", "rel.dev.(%)",
                "baseline(%)");
    for (const auto& baseline : lcdBaselines()) {
        const auto structure = aleator::io::readCif(poreGeometryDataFile(baseline.name + ".cif"));
        const std::vector<double> radii = radiiFor(structure.speciesSymbols);

        PoreAnalysisOptions options;
        options.asaSamplesPerAtom = 1; // LCD/PLD don't use MC sampling; keep this cheap
        options.volumeSamplesTotal = 1;

        const PoreGeometry result =
            analyzePoreGeometry(structure.particles, structure.lattice, radii, options);

        const double relDeviation =
            std::abs(result.largestCavityDiameterAngstrom - baseline.zeoPlusPlusLcdAngstrom) /
            baseline.zeoPlusPlusLcdAngstrom;
        const double allowedMax = baseline.baselineRelativeDeviation + kAllowedDeviationGrowth;

        std::printf("%10s %14.5f %14.5f %14.2f %14.2f\n", baseline.name.c_str(),
                    result.largestCavityDiameterAngstrom, baseline.zeoPlusPlusLcdAngstrom,
                    100.0 * relDeviation, 100.0 * baseline.baselineRelativeDeviation);

        INFO(baseline.name << ": LCD relative deviation=" << 100.0 * relDeviation
                           << "%, baseline=" << 100.0 * baseline.baselineRelativeDeviation
                           << "%, allowed max=" << 100.0 * allowedMax << "%");
        CHECK(relDeviation < allowedMax);
    }
}

TEST_CASE("PTY's PLD/ASA/AV deviation from real Zeo++ 0.4.7 has not grown beyond its recorded "
          "baseline",
          "[known-deviation][geometry][pore][slow]") {
    const auto structure = aleator::io::readCif(poreGeometryDataFile("PTY.cif"));
    const std::vector<double> radii = radiiFor(structure.speciesSymbols);

    PoreAnalysisOptions options;
    options.probeRadiusAngstrom = 1.86;

    const PoreGeometry result =
        analyzePoreGeometry(structure.particles, structure.lattice, radii, options);

    constexpr double kZeoPlusPlusPld = 3.92186;
    constexpr double kZeoPlusPlusTotalArea = 45.0748; // ASA + NASA, NASA=0
    constexpr double kZeoPlusPlusTotalVolume = 9.92047; // AV + NAV, NAV=0

    const double totalArea =
        result.accessibleSurfaceAreaAngstromSq + result.inaccessibleSurfaceAreaAngstromSq;
    const double totalVolume =
        result.accessibleVolumeAngstromCubed + result.inaccessibleVolumeAngstromCubed;

    const double pldRelDeviation =
        std::abs(result.poreLimitingDiameterAngstrom - kZeoPlusPlusPld) / kZeoPlusPlusPld;
    const double areaRelDeviation = std::abs(totalArea - kZeoPlusPlusTotalArea) / kZeoPlusPlusTotalArea;
    const double volumeRelDeviation =
        std::abs(totalVolume - kZeoPlusPlusTotalVolume) / kZeoPlusPlusTotalVolume;

    std::printf("PTY: PLD computed=%.5f Zeo++=%.5f (%.2f%%) | area computed=%.3f Zeo++=%.3f (%.2f%%) "
                "| volume computed=%.4f Zeo++=%.4f (%.2f%%)\n",
                result.poreLimitingDiameterAngstrom, kZeoPlusPlusPld, 100.0 * pldRelDeviation, totalArea,
                kZeoPlusPlusTotalArea, 100.0 * areaRelDeviation, totalVolume, kZeoPlusPlusTotalVolume,
                100.0 * volumeRelDeviation);

    // Baselines recorded from a real run of this exact test, at the same
    // PoreAnalysisOptions defaults used here (asaSamplesPerAtom=2000,
    // volumeSamplesTotal=50000, sampleSeed=0 -- deterministic, CLAUDE.md
    // invariant #5, so this reproduces bit-for-bit): PLD 1.81%, area
    // 11.98%, volume 8.78%. See known_deviation_baseline.md. Margins sized
    // the same way as the LCD test above.
    INFO("PTY PLD relative deviation=" << 100.0 * pldRelDeviation << "%");
    CHECK(pldRelDeviation < 0.018 + kAllowedDeviationGrowth);
    INFO("PTY total area relative deviation=" << 100.0 * areaRelDeviation << "%");
    CHECK(areaRelDeviation < 0.12 + kAllowedDeviationGrowth);
    INFO("PTY total volume relative deviation=" << 100.0 * volumeRelDeviation << "%");
    CHECK(volumeRelDeviation < 0.088 + kAllowedDeviationGrowth);
}
