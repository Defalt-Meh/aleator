// Validation anchor for the mixture-GCMC milestone's item 1, the gating
// check that must pass before anything else in the milestone proceeds
// ("Do not start (4) until (1) and (2) pass"):
//
//   "DEGENERATE CASE. A two-component simulation where both species have
//   identical parameters must reproduce the single-component isotherm at
//   the summed fugacity, and must give selectivity 1.0 within noise. This
//   catches most prefactor errors immediately and costs nothing to run."
//
// Two independent claims, both checked here against a real MonteCarloEngine
// run (not just the acceptance-ratio algebra already covered by
// tests/validation/test_gcmc_detailed_balance.cc):
//   1. <N_A> + <N_B> (two identical species, each at fugacityTotal/2) agrees
//      with <N_ref> (one species, at fugacityTotal) within combined
//      statistical error.
//   2. Selectivity S = (x_A/x_B)/(y_A/y_B) = N_A/N_B (since y_A=y_B=0.5 by
//      construction: fugacityTotal split evenly) is 1.0 within noise --
//      identity-swap moves are what makes this equilibrate at a practical
//      rate (pure insertion/deletion alone would take much longer to mix
//      the two labels evenly, per the milestone's own note).
//
// A small synthetic framework (not a real CIF) is used deliberately, same
// convention as test_gcmc_detailed_balance.cc/test_gcmc_cell_list_bit_identical.cc:
// this test is about mixture bookkeeping/acceptance correctness, not about
// reproducing a specific real adsorption isotherm (that's a separate,
// later item in the milestone, gated on this one passing first).

#include <catch2/catch_test_macros.hpp>

#include <array>
#include <cmath>
#include <cstdio>
#include <memory>
#include <random>
#include <vector>

#include "core/geometry/lattice.hpp"
#include "core/math/counter_based_rng.hpp"
#include "core/math/particle_data.hpp"
#include "engines/monte_carlo/molecule_species.hpp"
#include "engines/monte_carlo/monte_carlo_engine.hpp"
#include "forcefield/pairwise/lennard_jones.hpp"

using aleator::core::Lattice;
using aleator::core::ParticleData;
using aleator::core::Philox4x32Rng;
using aleator::engines::MoleculeSpecies;
using aleator::engines::MonteCarloEngine;
using aleator::forcefield::LennardJones;
using aleator::forcefield::LennardJonesParameters;
using aleator::forcefield::LennardJonesTruncation;

namespace {

constexpr double kBoxLengthAngstrom = 25.0;
constexpr double kCutoffAngstrom = 10.0;
constexpr std::size_t kFrameworkAtomCount = 120;
constexpr double kTemperatureKelvin = 298.0;

Lattice cubicLattice() {
    std::array<std::array<double, 3>, 3> m{};
    m[0] = {kBoxLengthAngstrom, 0.0, 0.0};
    m[1] = {0.0, kBoxLengthAngstrom, 0.0};
    m[2] = {0.0, 0.0, kBoxLengthAngstrom};
    return Lattice(m);
}

ParticleData randomFramework(const Lattice& lattice, unsigned seed) {
    std::mt19937 rng(seed);
    std::uniform_real_distribution<double> unit(0.0, 1.0);
    ParticleData framework;
    framework.resize(kFrameworkAtomCount);
    for (std::size_t i = 0; i < kFrameworkAtomCount; ++i) {
        const auto cart = lattice.fractionalToCartesian({unit(rng), unit(rng), unit(rng)});
        framework.x[i] = cart[0];
        framework.y[i] = cart[1];
        framework.z[i] = cart[2];
        framework.mass[i] = 40.0;
        framework.charge[i] = 0.0;
        framework.species[i] = 0;
    }
    return framework;
}

struct BlockStats {
    double mean = 0.0;
    double standardError = 0.0;
};

/// Runs `equilibrationSteps` (off any accounting) then `productionSteps`
/// split into `blocks` blocks, sampling `sample(engine)` after every
/// production step and block-averaging -- the same estimator/error-bar
/// methodology test_gcmc_ch4_irmof1_isotherm.cc and
/// test_gcmc_co2_irmof1_isotherm.cc already use.
template <typename SampleFn>
BlockStats runBlocked(MonteCarloEngine& engine, std::size_t equilibrationSteps,
                       std::size_t productionSteps, int blocks, SampleFn sample) {
    engine.run(equilibrationSteps);
    const std::size_t blockSteps = productionSteps / static_cast<std::size_t>(blocks);
    std::vector<double> blockMeans;
    blockMeans.reserve(static_cast<std::size_t>(blocks));
    for (int b = 0; b < blocks; ++b) {
        double sum = 0.0;
        for (std::size_t step = 0; step < blockSteps; ++step) {
            engine.run(1);
            sum += sample(engine);
        }
        blockMeans.push_back(sum / static_cast<double>(blockSteps));
    }
    double mean = 0.0;
    for (double v : blockMeans) {
        mean += v;
    }
    mean /= static_cast<double>(blocks);
    double variance = 0.0;
    for (double v : blockMeans) {
        variance += (v - mean) * (v - mean);
    }
    variance /= static_cast<double>(blocks - 1);
    return {mean, std::sqrt(variance / static_cast<double>(blocks))};
}

} // namespace

TEST_CASE("Mixture GCMC degenerate case: two identical species reproduce the single-species "
          "isotherm at the summed fugacity, and give selectivity 1.0",
          "[validation][montecarlo][mixture]") {
    const Lattice lattice = cubicLattice();
    const auto framework = randomFramework(lattice, 20260812);

    // Framework species 0, guest LJ species 1 -- both species A and B use
    // the SAME LJ species index (identical parameters is the whole point of
    // this test), so there is exactly one guest LJ entry.
    const std::vector<LennardJonesParameters> ljParameters{
        {50.0, 3.0},   // framework, species 0
        {120.0, 3.5},  // guest, species 1 (both A and B use this)
    };
    auto forceField = std::make_shared<const LennardJones>(ljParameters, kCutoffAngstrom,
                                                             LennardJonesTruncation::Shifted);

    MoleculeSpecies guest;
    guest.sites.push_back({0.0, 0.0, 0.0, /*ljSpecies=*/1, /*charge=*/0.0, /*mass=*/16.0});

    // Calibrated empirically (this file's own value, per CLAUDE.md
    // invariant #1 -- a first try at 3.5e-3 gave only <N> ~ 0.16, too low
    // for insertion/deletion/swap to have real material to act on; 0.3
    // gives <N> ~ 12-13 total, ~6-7 per species): a real, moderate
    // occupancy at this system's size.
    constexpr double kFugacityTotalInternal = 0.3;
    constexpr double kPascalToInternal = MonteCarloEngine::kPascalToInternal;
    const double fugacityTotalPascal = kFugacityTotalInternal / kPascalToInternal;

    constexpr std::size_t kEquilibrationSteps = 40'000;
    constexpr std::size_t kProductionSteps = 200'000;
    constexpr int kBlocks = 10;

    // Reference: one species, full fugacity.
    double refMean;
    double refError;
    {
        auto rng = std::make_unique<Philox4x32Rng>();
        rng->seed(/*seed=*/7001, /*streamIndex=*/0);
        MonteCarloEngine engine(framework, lattice, forceField, std::move(rng), kTemperatureKelvin,
                                 guest, fugacityTotalPascal);
        const auto stats = runBlocked(engine, kEquilibrationSteps, kProductionSteps, kBlocks,
                                       [](const MonteCarloEngine& e) {
                                           return static_cast<double>(e.moleculeCount());
                                       });
        refMean = stats.mean;
        refError = stats.standardError;
    }

    // Mixture: species A and B, identical parameters, each at half the
    // total fugacity.
    double mixtureTotalMean;
    double mixtureTotalError;
    double nAMean;
    double nAError;
    double nBMean;
    double nBError;
    {
        auto rng = std::make_unique<Philox4x32Rng>();
        rng->seed(/*seed=*/7002, /*streamIndex=*/0);
        MonteCarloEngine engine(framework, lattice, forceField, std::move(rng), kTemperatureKelvin,
                                 std::vector<MoleculeSpecies>{guest, guest},
                                 std::vector<double>{fugacityTotalPascal / 2.0,
                                                      fugacityTotalPascal / 2.0});
        // Sampled directly here (not via runBlocked(), which only samples
        // one quantity per call) so total/A/B counts come from the exact
        // same trajectory rather than three separately-seeded runs.
        engine.run(kEquilibrationSteps);
        const std::size_t blockSteps = kProductionSteps / static_cast<std::size_t>(kBlocks);
        std::vector<double> totalMeans, aMeans, bMeans;
        for (int b = 0; b < kBlocks; ++b) {
            double sumTotal = 0.0, sumA = 0.0, sumB = 0.0;
            for (std::size_t step = 0; step < blockSteps; ++step) {
                engine.run(1);
                const double nA = static_cast<double>(engine.moleculeCountOfSpecies(0));
                const double nB = static_cast<double>(engine.moleculeCountOfSpecies(1));
                sumTotal += nA + nB;
                sumA += nA;
                sumB += nB;
            }
            totalMeans.push_back(sumTotal / static_cast<double>(blockSteps));
            aMeans.push_back(sumA / static_cast<double>(blockSteps));
            bMeans.push_back(sumB / static_cast<double>(blockSteps));
        }
        auto meanAndError = [&](const std::vector<double>& v) -> BlockStats {
            double m = 0.0;
            for (double x : v) m += x;
            m /= static_cast<double>(v.size());
            double var = 0.0;
            for (double x : v) var += (x - m) * (x - m);
            var /= static_cast<double>(v.size() - 1);
            return {m, std::sqrt(var / static_cast<double>(v.size()))};
        };
        const auto totalStatsReal = meanAndError(totalMeans);
        const auto aStats = meanAndError(aMeans);
        const auto bStats = meanAndError(bMeans);
        mixtureTotalMean = totalStatsReal.mean;
        mixtureTotalError = totalStatsReal.standardError;
        nAMean = aStats.mean;
        nAError = aStats.standardError;
        nBMean = bStats.mean;
        nBError = bStats.standardError;
    }

    const double combinedTotalError = std::sqrt(refError * refError + mixtureTotalError * mixtureTotalError);
    const double totalSigmas =
        combinedTotalError > 0.0 ? std::abs(mixtureTotalMean - refMean) / combinedTotalError : 0.0;

    // Selectivity S = (x_A/x_B) / (y_A/y_B); y_A=y_B=0.5 by construction
    // (fugacity split evenly), so S = x_A/x_B = N_A/N_B directly.
    const double selectivity = nAMean / nBMean;
    // Error propagation for a ratio of two independent means:
    // sigma_S/S = sqrt((sigma_A/N_A)^2 + (sigma_B/N_B)^2).
    const double selectivityRelError =
        std::sqrt((nAError / nAMean) * (nAError / nAMean) + (nBError / nBMean) * (nBError / nBMean));
    const double selectivityError = selectivity * selectivityRelError;
    const double selectivitySigmas =
        selectivityError > 0.0 ? std::abs(selectivity - 1.0) / selectivityError : 0.0;

    std::printf("reference (1 species, full fugacity): <N> = %.3f +/- %.3f\n", refMean, refError);
    std::printf("mixture (2 identical species, split fugacity): <N_A+N_B> = %.3f +/- %.3f "
                "(%.2f sigma from reference)\n",
                mixtureTotalMean, mixtureTotalError, totalSigmas);
    std::printf("  <N_A> = %.3f +/- %.3f, <N_B> = %.3f +/- %.3f\n", nAMean, nAError, nBMean, nBError);
    std::printf("  selectivity S = N_A/N_B = %.4f +/- %.4f (%.2f sigma from 1.0)\n", selectivity,
                selectivityError, selectivitySigmas);

    INFO("refMean=" << refMean << " +/- " << refError << ", mixtureTotalMean=" << mixtureTotalMean
                     << " +/- " << mixtureTotalError << ", sigmas=" << totalSigmas);
    REQUIRE(totalSigmas < 3.0);

    INFO("selectivity=" << selectivity << " +/- " << selectivityError
                         << ", sigmas from 1.0=" << selectivitySigmas);
    REQUIRE(selectivitySigmas < 3.0);
}
