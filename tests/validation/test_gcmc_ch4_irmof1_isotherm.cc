// Validation anchor (CLAUDE.md milestone: "Reproduce a published CH4
// adsorption isotherm in IRMOF-1... Do not declare it passing unless the
// curve genuinely overlays within error.") and CLAUDE.md section 4 ("Handling
// a validated component that disagrees with its reference"): this file used
// to contain a single test comparing the full computed curve against pyIAST's
// published isotherm on a tolerance band wide enough (~20%, vs. an observed
// ~14.5% worst case) to swallow a genuine 10% regression introduced
// tomorrow. Per section 4 it has been split in two:
//
//   - THIS test: tight, self-contained, no external reference. It checks
//     what this codebase's own internals genuinely pin down -- that the
//     real IRMOF-1/methane system's low-pressure loading is consistent with
//     this same engine's own Widom-insertion Henry coefficient on the same
//     real structure and force field. A regression in the acceptance
//     criteria, the force-field call, or the fugacity conversion would very
//     likely break this, independent of whether the absolute loading
//     matches the literature.
//   - tests/known_deviation/test_gcmc_ch4_irmof1_known_deviation.cc: the
//     informational full-curve comparison against the published pyIAST
//     reference, which only fails if the gap widens beyond a checked-in
//     baseline (tests/validation/data/irmof1/known_deviation_baseline.md).
//
// Real structure (RASPA2's IRMOF-1.cif, symmetry-expanded by this
// codebase's own io/ CIF reader), real force field (RASPA2's GenericMOFs
// UFF-derived framework parameters + TraPPE-style CH4_sp3 united-atom
// adsorbate). See data/irmof1/PROVENANCE.md for exact sources.

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
    return std::filesystem::path(ALEATOR_TEST_DATA_DIR) / "data" / "irmof1" / name;
}

/// UFF-derived LJ parameters from RASPA2's "GenericMOFs" force field (see
/// data/irmof1/PROVENANCE.md) — epsilon in K, sigma in Angstrom, already in
/// this codebase's internal unit convention. Re-verified live against
/// RASPA2's force_field_mixing_rules.def and pseudo_atoms.def during the
/// CLAUDE.md hygiene session that split this test (not just trusted from an
/// earlier session's transcription).
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

} // namespace

TEST_CASE("GCMC-computed CH4/IRMOF-1 loading at low pressure is self-consistent with this "
          "engine's own Widom-insertion Henry coefficient",
          "[validation][montecarlo][slow]") {
    const auto structure = aleator::io::readCif(irmof1DataFile("IRMOF-1.cif"));
    const std::size_t frameworkCount = structure.particles.size();
    REQUIRE(frameworkCount > 0);

    std::vector<LennardJonesParameters> ljParameters;
    for (const auto& symbol : structure.speciesSymbols) {
        ljParameters.push_back(uffParameters(symbol));
    }
    const auto ch4SpeciesIndex = static_cast<std::uint32_t>(ljParameters.size());
    ljParameters.push_back({158.5, 3.72}); // CH4_sp3, RASPA2 GenericMOFs

    const double cutoff = 12.0; // < L_perp/2 (25.832/2 = 12.916 Ang for this cubic cell)
    auto forceField = std::make_shared<const LennardJones>(ljParameters, cutoff,
                                                             LennardJonesTruncation::Shifted);

    MoleculeSpecies ch4;
    ch4.sites.push_back({0.0, 0.0, 0.0, ch4SpeciesIndex, /*charge=*/0.0, /*mass=*/16.04246});

    const double temperatureKelvin = 298.0;
    const PengRobinson eos(kMethane);
    const double volume = std::abs(structure.lattice.volume());

    // Widom insertion against the (unpopulated) real framework: cheap, done
    // once, no external reference involved.
    auto widomRng = std::make_unique<Philox4x32Rng>();
    widomRng->seed(/*seed=*/2026, /*streamIndex=*/999);
    auto widomParticles = structure.particles;
    MonteCarloEngine widomEngine(widomParticles, structure.lattice, forceField,
                                  std::move(widomRng), temperatureKelvin, ch4,
                                  /*fugacityPascal=*/1.0);
    const double henryCoefficient = widomEngine.widomInsertionHenryCoefficient(500'000);
    std::printf("Real IRMOF-1 Widom insertion Henry coefficient K_H = %.6e (1/K)\n",
                henryCoefficient);
    REQUIRE(henryCoefficient > 0.0);

    // A real GCMC production run at the same low pressure this codebase's
    // known-deviation test's lowest point uses (0.1 bar) -- deep in the
    // dilute/Henry regime for this system (see the known-deviation test's
    // printed <N> ~0.2-0.3 out of a system that reaches <N> ~23 by 10 bar).
    const double pressurePascal = 0.1 * 1.0e5;
    const double fugacityPascal = eos.fugacityPascal(temperatureKelvin, pressurePascal);
    const double fugacityInternal = fugacityPascal * MonteCarloEngine::kPascalToInternal;

    constexpr std::size_t kEquilibrationSteps = 30'000;
    constexpr std::size_t kProductionSteps = 90'000;
    constexpr int kBlocks = 10;
    constexpr std::size_t kBlockSteps = kProductionSteps / kBlocks;

    auto rng = std::make_unique<Philox4x32Rng>();
    rng->seed(/*seed=*/2026, /*streamIndex=*/0);
    auto particles = structure.particles;
    MonteCarloEngine engine(particles, structure.lattice, forceField, std::move(rng),
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

    const double loading = meanN / volume;
    const double loadingStandardError = standardErrorN / volume;
    const double impliedHenryCoefficient = loading / fugacityInternal;
    const double relError =
        std::abs(impliedHenryCoefficient - henryCoefficient) / henryCoefficient;
    const double relStandardError = loadingStandardError / fugacityInternal / henryCoefficient;

    std::printf(
        "P=0.1 bar: <N> = %.4f +/- %.4f, GCMC-implied K_H = %.6e +/- %.6e (1/K) "
        "[%.2f%% rel. std. error]\nWidom K_H          = %.6e (1/K)\n"
        "relative difference = %.3f%%\n",
        meanN, standardErrorN, impliedHenryCoefficient, loadingStandardError / fugacityInternal,
        100.0 * relStandardError, henryCoefficient, 100.0 * relError);

    // Tight: this is a self-consistency check against this codebase's own
    // Henry coefficient on the SAME real structure, not against any
    // external published number, so there is no systematic-offset excuse
    // available here the way there is for the known-deviation test.
    // Measured during development: relative difference 3.42%, block-
    // averaged relative standard error 1.49% (~2.3 standard errors) --
    // i.e. the sampling machinery on the real IRMOF-1 system already
    // agrees with its own Henry-limit prediction far more closely than it
    // agrees with the external pyIAST reference (12-15%, known_deviation
    // test), which is itself evidence the ~12-15% gap is not a sampling
    // bug. 12% leaves a ~3.5x margin over the measured 3.42% (room for
    // legitimate cross-compiler/cross-platform statistical noise in a
    // chaotic MC trajectory -- CLAUDE.md invariant #5's bit-identical
    // reproducibility is scoped to a fixed build, not guaranteed across
    // compilers) while still being tight enough that a real 10-point
    // regression in the sampling machinery would fail it (CLAUDE.md
    // definition of done).
    constexpr double kMaxRelativeError = 0.12;
    REQUIRE(relError < kMaxRelativeError);
}
