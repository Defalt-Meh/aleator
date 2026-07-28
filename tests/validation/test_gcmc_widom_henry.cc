// Validation anchor (CLAUDE.md milestone: "Widom insertion Henry
// coefficient must match the low-pressure limiting slope of a computed
// isotherm (self-consistency, 2%)"). Self-contained synthetic system (a
// small fixed LJ "framework" plus a single-site LJ "adsorbate" species) —
// deliberately NOT the real IRMOF-1/methane system (that's
// test_gcmc_ch4_irmof1_isotherm.cc); this test isolates the GCMC engine's
// own internal consistency (does its equilibrium loading at low fugacity
// actually reproduce what Widom insertion predicts?) from any question of
// matching an external published curve.
//
// Statistics: several independent replicas (different RNG streams) rather
// than one long run, so the reported uncertainty is an honest
// standard-error-of-the-mean across genuinely independent trajectories,
// not a guess at the single-run autocorrelation time.

#include <catch2/catch_test_macros.hpp>

#include <array>
#include <cmath>
#include <cstdio>
#include <memory>
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

Lattice cubicLattice(double length) {
    std::array<std::array<double, 3>, 3> m{};
    m[0] = {length, 0.0, 0.0};
    m[1] = {0.0, length, 0.0};
    m[2] = {0.0, 0.0, length};
    return Lattice(m);
}

ParticleData buildFramework(double length) {
    ParticleData framework;
    const std::array<double, 2> coords{length * 0.25, length * 0.75};
    for (double x : coords) {
        for (double y : coords) {
            for (double z : coords) {
                framework.push_back(x, y, z, 0.0, 0.0, 0.0, /*mass=*/65.0, /*charge=*/0.0,
                                     /*species=*/0);
            }
        }
    }
    return framework;
}

} // namespace

TEST_CASE("GCMC low-pressure isotherm slope matches Widom-insertion Henry coefficient",
          "[validation][montecarlo]") {
    const double boxLength = 16.0;
    const Lattice lattice = cubicLattice(boxLength);
    const ParticleData framework = buildFramework(boxLength);

    const std::vector<LennardJonesParameters> species{
        {100.0, 3.0}, // 0: framework
        {100.0, 3.5}, // 1: single-site adsorbate
    };
    const double cutoff = 6.0; // < L_perp/2 = 8.0
    auto forceField = std::make_shared<const LennardJones>(species, cutoff,
                                                             LennardJonesTruncation::Shifted);

    MoleculeSpecies adsorbate;
    adsorbate.sites.push_back({0.0, 0.0, 0.0, /*ljSpecies=*/1, /*charge=*/0.0, /*mass=*/16.0});

    const double temperature = 150.0;

    // Widom insertion against the (unpopulated) framework: cheap, done once.
    auto widomRng = std::make_unique<Philox4x32Rng>();
    widomRng->seed(/*seed=*/1, /*streamIndex=*/0);
    MonteCarloEngine widomEngine(framework, lattice, forceField, std::move(widomRng), temperature,
                                  adsorbate, /*fugacityPascal=*/1.0);
    const double henryCoefficient = widomEngine.widomInsertionHenryCoefficient(2'000'000);
    std::printf("Widom-insertion Henry coefficient K_H = %.6e (1/K)\n", henryCoefficient);
    REQUIRE(henryCoefficient > 0.0);

    // Pick a fugacity low enough that predicted mean loading (K_H * f * V)
    // stays well under 1 molecule -- the low-pressure / dilute regime the
    // Henry coefficient is defined for.
    const double volume = std::abs(lattice.volume());
    const double targetMeanN = 0.3;
    const double fugacityInternal = targetMeanN / (henryCoefficient * volume);
    const double fugacityPascal = fugacityInternal / MonteCarloEngine::kPascalToInternal;
    std::printf("Target <N>=%.3f -> fugacityInternal=%.6e K/Ang^3, fugacityPascal=%.6e Pa\n",
                targetMeanN, fugacityInternal, fugacityPascal);

    constexpr int kReplicas = 6;
    constexpr std::size_t kEquilibrationSteps = 40'000;
    constexpr std::size_t kProductionSteps = 200'000;
    constexpr std::size_t kSampleStride = 20;

    std::vector<double> replicaLoadings;
    for (int replica = 0; replica < kReplicas; ++replica) {
        auto rng = std::make_unique<Philox4x32Rng>();
        rng->seed(/*seed=*/1000, /*streamIndex=*/static_cast<std::uint64_t>(replica));
        MonteCarloEngine engine(framework, lattice, forceField, std::move(rng), temperature,
                                 adsorbate, fugacityPascal);

        engine.run(kEquilibrationSteps);

        double sumN = 0.0;
        std::size_t samples = 0;
        for (std::size_t step = 0; step < kProductionSteps; step += kSampleStride) {
            engine.run(kSampleStride);
            sumN += static_cast<double>(engine.moleculeCount());
            ++samples;
        }
        const double meanN = sumN / static_cast<double>(samples);
        replicaLoadings.push_back(meanN / volume); // molecules / Angstrom^3
        std::printf("  replica %d: <N> = %.4f, loading = %.6e molecules/Ang^3\n", replica, meanN,
                    replicaLoadings.back());
    }

    double meanLoading = 0.0;
    for (double l : replicaLoadings) {
        meanLoading += l;
    }
    meanLoading /= static_cast<double>(kReplicas);

    double variance = 0.0;
    for (double l : replicaLoadings) {
        variance += (l - meanLoading) * (l - meanLoading);
    }
    variance /= static_cast<double>(kReplicas - 1);
    const double standardErrorOfMean = std::sqrt(variance / static_cast<double>(kReplicas));

    const double impliedHenryCoefficient = meanLoading / fugacityInternal;
    const double relError =
        std::abs(impliedHenryCoefficient - henryCoefficient) / henryCoefficient;
    const double relStandardError = standardErrorOfMean / fugacityInternal / henryCoefficient;

    std::printf(
        "GCMC-implied K_H = %.6e +/- %.6e (1/K) [%.2f%% rel. std. error]\n"
        "Widom K_H         = %.6e (1/K)\n"
        "relative difference = %.3f%%\n",
        impliedHenryCoefficient, standardErrorOfMean / fugacityInternal,
        100.0 * relStandardError, henryCoefficient, 100.0 * relError);

    REQUIRE(relError < 0.02);
}
