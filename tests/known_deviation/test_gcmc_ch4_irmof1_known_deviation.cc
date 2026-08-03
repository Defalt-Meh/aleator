// Informational known-deviation test (CLAUDE.md section 4: "Handling a
// validated component that disagrees with its reference"). This is the
// full four-point CH4/IRMOF-1 curve compared against the published pyIAST
// reference isotherm (Simon/Smit/Haranczyk, Computer Physics
// Communications 2016) -- moved out of tests/validation/ on purpose. This
// codebase's own sampling machinery is independently, tightly verified
// elsewhere (exact detailed balance:
// tests/validation/test_gcmc_detailed_balance.cc; Widom self-consistency
// on this exact real system: tests/validation/test_gcmc_ch4_irmof1_isotherm.cc)
// -- this test instead tracks whether the REMAINING, still-unexplained gap
// against the external reference is growing. It is deliberately NOT a tight
// pass/fail against pyIAST's numbers: CLAUDE.md section 4 is explicit that
// widening a tolerance to swallow an unexplained disagreement is a defect
// disguised as a passing test, so this test instead compares against a
// checked-in baseline (tests/validation/data/irmof1/
// known_deviation_baseline.md) and fails only if the deviation has grown.
//
// Investigation performed and what was ruled in/out: see
// tests/validation/data/irmof1/known_deviation_baseline.md. Short version:
// tail-correction convention, mixing rule, and LJ parameters were all
// checked against RASPA2's actual force-field files (re-fetched live, not
// assumed) and match exactly; pore blocking was ruled out by reading this
// codebase's own insertion code (uniform full-cell sampling, no
// accessibility filtering); a real sensitivity check ruled out the
// within-cell cutoff margin (12.0 vs. 12.8 Ang at P=1 bar moved <N> by
// only 0.25 standard errors) but NOT a larger cutoff enabled by a bigger
// simulation cell, which was not attempted (too slow at this codebase's
// current O(N)-scan GCMC engine, CLAUDE.md section 0 defect 4); pyIAST's
// own (undocumented) simulation parameters remain unverifiable from
// available sources. This is an owned, unresolved defect (CLAUDE.md
// section 0), not closed here.

#include <catch2/catch_test_macros.hpp>

#include <cmath>
#include <cstdio>
#include <filesystem>
#include <map>
#include <memory>
#include <stdexcept>
#include <string>
#include <vector>

#include "core/math/counter_based_rng.hpp"
#include "engines/monte_carlo/molecule_species.hpp"
#include "engines/monte_carlo/monte_carlo_engine.hpp"
#include "engines/monte_carlo/peng_robinson.hpp"
#include "forcefield/pairwise/lennard_jones.hpp"
#include "io/structure_io.hpp"

using aleator::core::Philox4x32Rng;
using aleator::engines::kMethane;
using aleator::engines::MoleculeSpecies;
using aleator::engines::MonteCarloEngine;
using aleator::engines::PengRobinson;
using aleator::forcefield::LennardJones;
using aleator::forcefield::LennardJonesParameters;
using aleator::forcefield::LennardJonesTruncation;

namespace {

std::filesystem::path irmof1DataFile(const std::string& name) {
    return std::filesystem::path(ALEATOR_KNOWN_DEVIATION_DATA_DIR) / "data" / "irmof1" / name;
}

LennardJonesParameters uffParameters(const std::string& element) {
    static const std::map<std::string, LennardJonesParameters> table{
        {"Zn", {62.3992, 2.46155}},
        {"O", {48.1581, 3.03315}},
        {"C", {47.8562, 3.47299}},
        {"H", {7.64893, 2.84642}},
    };
    const auto it = table.find(element);
    if (it == table.end()) {
        throw std::runtime_error("uffParameters: no GenericMOFs LJ parameters for element '" +
                                  element + "'");
    }
    return it->second;
}

struct ReferencePoint {
    double pressureBar;
    double loadingMmolPerG;
    /// Recorded from a real run of this test (see known_deviation_baseline.md
    /// for the full provenance/investigation writeup) -- the maximum
    /// relative deviation this test tolerates before failing is derived
    /// from this, not from the reference value alone.
    double baselineRelativeDeviation;
};

/// Transcribed directly from data/irmof1/methane_irmof1_isotherm_298K.csv
/// (see PROVENANCE.md); baselineRelativeDeviation transcribed from
/// known_deviation_baseline.md.
const std::vector<ReferencePoint>& referenceIsotherm() {
    static const std::vector<ReferencePoint> data{
        {0.1, 0.04301247025401385, 0.1258},
        {1.0, 0.4384412654666955, 0.1196},
        {5.0, 2.2082018401014736, 0.1368},
        {10.0, 4.389245423072791, 0.1445},
    };
    return data;
}

constexpr double kAmuToKg = 1.66053906660e-27; // CODATA 2018
constexpr double kAvogadro = 6.02214076e23;    // SI 2019 exact

// Additive margin (in relative-deviation percentage points) allowed on top
// of each point's recorded baseline before this test fails. Sized to
// absorb re-run-to-re-run Monte Carlo noise (the block-averaged standard
// error on loading measured here is ~1.5-2% relative -- see printed
// output), with real room left over so a genuine several-point widening
// still fails. Not sized to also swallow a regression -- CLAUDE.md
// section 4: "never widen a tolerance to accommodate a disagreement you
// have not explained."
constexpr double kAllowedDeviationGrowth = 0.05;

} // namespace

TEST_CASE("GCMC-computed CH4/IRMOF-1 isotherm's deviation from the published pyIAST reference "
          "has not grown beyond its recorded baseline",
          "[known-deviation][montecarlo][slow]") {
    const auto structure = aleator::io::readCif(irmof1DataFile("IRMOF-1.cif"));
    const std::size_t frameworkCount = structure.particles.size();
    REQUIRE(frameworkCount > 0);

    std::vector<LennardJonesParameters> ljParameters;
    for (const auto& symbol : structure.speciesSymbols) {
        ljParameters.push_back(uffParameters(symbol));
    }
    const auto ch4SpeciesIndex = static_cast<std::uint32_t>(ljParameters.size());
    ljParameters.push_back({158.5, 3.72}); // CH4_sp3, RASPA2 GenericMOFs

    double frameworkMassAmu = 0.0;
    for (double m : structure.particles.mass) {
        frameworkMassAmu += m;
    }
    const double frameworkMassKg = frameworkMassAmu * kAmuToKg;
    const double volumeAngstrom3 = std::abs(structure.lattice.volume());

    std::printf("IRMOF-1: %zu framework atoms, %zu species, cell volume = %.2f Ang^3, "
                "framework mass = %.3f amu\n",
                frameworkCount, structure.speciesSymbols.size(), volumeAngstrom3, frameworkMassAmu);

    const double cutoff = 12.0; // < L_perp/2 (25.832/2 = 12.916 Ang for this cubic cell)
    auto forceField = std::make_shared<const LennardJones>(ljParameters, cutoff,
                                                             LennardJonesTruncation::Shifted);

    MoleculeSpecies ch4;
    ch4.sites.push_back({0.0, 0.0, 0.0, ch4SpeciesIndex, /*charge=*/0.0, /*mass=*/16.04246});

    const double temperatureKelvin = 298.0;
    const PengRobinson eos(kMethane);

    constexpr std::size_t kEquilibrationSteps = 30'000;
    constexpr std::size_t kProductionSteps = 90'000;
    constexpr int kBlocks = 10;
    constexpr std::size_t kBlockSteps = kProductionSteps / kBlocks;

    std::vector<double> computedLoadings;
    std::vector<double> computedErrors;

    for (std::size_t pointIndex = 0; pointIndex < referenceIsotherm().size(); ++pointIndex) {
        const auto& ref = referenceIsotherm()[pointIndex];
        const double pressurePascal = ref.pressureBar * 1.0e5;
        const double fugacityPascal = eos.fugacityPascal(temperatureKelvin, pressurePascal);

        auto rng = std::make_unique<Philox4x32Rng>();
        rng->seed(/*seed=*/2026, /*streamIndex=*/static_cast<std::uint64_t>(pointIndex));
        MonteCarloEngine engine(structure.particles, structure.lattice, forceField, std::move(rng),
                                 temperatureKelvin, ch4, fugacityPascal);

        engine.run(kEquilibrationSteps);

        std::vector<double> blockMeans;
        blockMeans.reserve(kBlocks);
        for (int block = 0; block < kBlocks; ++block) {
            double sumN = 0.0;
            for (std::size_t step = 0; step < kBlockSteps; ++step) {
                engine.run(1);
                sumN += static_cast<double>(engine.moleculeCount());
            }
            blockMeans.push_back(sumN / static_cast<double>(kBlockSteps));
        }

        double meanN = 0.0;
        for (double b : blockMeans) {
            meanN += b;
        }
        meanN /= static_cast<double>(kBlocks);

        double blockVariance = 0.0;
        for (double b : blockMeans) {
            blockVariance += (b - meanN) * (b - meanN);
        }
        blockVariance /= static_cast<double>(kBlocks - 1);
        const double standardErrorN = std::sqrt(blockVariance / static_cast<double>(kBlocks));

        // mol/kg and mmol/g are numerically identical, so no extra factor
        // belongs here -- see known_deviation_baseline.md for the earlier
        // *1000.0 unit-conversion bug this exact computation caught.
        const double loadingMmolPerG = (meanN / kAvogadro) / frameworkMassKg;
        const double loadingErrorMmolPerG = (standardErrorN / kAvogadro) / frameworkMassKg;

        computedLoadings.push_back(loadingMmolPerG);
        computedErrors.push_back(loadingErrorMmolPerG);

        std::printf(
            "P=%7.3f bar  fugacity=%10.1f Pa  <N>=%7.3f +/- %6.3f  "
            "computed=%8.4f +/- %7.4f mmol/g  reference(pyIAST)=%8.4f mmol/g\n",
            ref.pressureBar, fugacityPascal, meanN, standardErrorN, loadingMmolPerG,
            loadingErrorMmolPerG, ref.loadingMmolPerG);
    }

    std::printf("\n--- Summary: computed vs. published (pyIAST GCMC, Simon/Smit/Haranczyk 2016) "
                "---\n");
    std::printf("%10s %18s %18s %14s %14s\n", "P (bar)", "computed (mmol/g)",
                "reference (mmol/g)", "rel. dev. (%)", "baseline (%)");
    for (std::size_t i = 0; i < referenceIsotherm().size(); ++i) {
        const double relDeviation =
            std::abs(computedLoadings[i] - referenceIsotherm()[i].loadingMmolPerG) /
            referenceIsotherm()[i].loadingMmolPerG;
        std::printf("%10.3f %18.4f %18.4f %14.2f %14.2f\n", referenceIsotherm()[i].pressureBar,
                    computedLoadings[i], referenceIsotherm()[i].loadingMmolPerG,
                    100.0 * relDeviation, 100.0 * referenceIsotherm()[i].baselineRelativeDeviation);
    }

    // Informational, not tight: fails only if the deviation from the
    // published reference has genuinely grown beyond its recorded
    // baseline (plus noise margin), never because the absolute gap itself
    // is nonzero -- that gap is the open, unresolved defect this test
    // exists to keep watch on, not to gate CI on directly (CLAUDE.md
    // section 4).
    for (std::size_t i = 0; i < referenceIsotherm().size(); ++i) {
        const double relDeviation =
            std::abs(computedLoadings[i] - referenceIsotherm()[i].loadingMmolPerG) /
            referenceIsotherm()[i].loadingMmolPerG;
        const double allowedMax =
            referenceIsotherm()[i].baselineRelativeDeviation + kAllowedDeviationGrowth;
        INFO("P=" << referenceIsotherm()[i].pressureBar
                   << " bar: relative deviation=" << 100.0 * relDeviation
                   << "%, baseline=" << 100.0 * referenceIsotherm()[i].baselineRelativeDeviation
                   << "%, allowed max=" << 100.0 * allowedMax << "%");
        CHECK(relDeviation < allowedMax);
    }
}
