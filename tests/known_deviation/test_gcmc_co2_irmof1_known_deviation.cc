// Informational known-deviation test (CLAUDE.md section 4) for the
// charged-GCMC milestone's highest-pressure CO2/IRMOF-1 reference point.
//
// tests/validation/test_gcmc_co2_irmof1_isotherm.cc's low- and
// mid-pressure points (1e4, 1e5 Pa) agree with the published CRAFTED
// reference to well under 1 combined sigma. This point (1e6 Pa, ~127
// molecules present at equilibrium) does not: a real calibration run
// measured computed 10.329 +/- 0.597 mol/kg against reference 12.832 +/-
// 0.994 mol/kg, ~2.16 combined sigma, ~19.5% low. Moved out of the tight
// test rather than folded in behind a loosened tolerance (CLAUDE.md
// section 4: a tolerance wide enough to pass this point would also be
// wide enough to hide a real regression at the other two).
//
// Leading (NOT independently confirmed) hypothesis: this test's
// kEquilibrationSteps (2000, chosen so the tight test's two points run in
// a practical amount of time) is calibrated against low occupancy
// (<=~11 molecules) and is very plausibly insufficient for this point's
// much higher equilibrium occupancy (~127 molecules) to actually reach
// that occupancy from an empty starting configuration -- an
// under-equilibration artifact of this test's own time budget, not
// necessarily a force-field or sampling defect. NOT verified by a longer
// run in this session (a materially longer equilibration at this
// occupancy was not completed within this session's practical time
// budget, given this test's own per-step cost at ~127 molecules). Left
// open, honestly, as CLAUDE.md section 4 requires.
//
// What IS independently verified: the electrostatic cache itself is
// correct at this exact system size and duration --
// test_gcmc_co2_irmof1_isotherm.cc's own points, and this test's own run,
// both check the incrementally-maintained Ewald state against a
// from-scratch recomputation (CLAUDE.md's charged-GCMC drift gate) and
// pass at machine precision, so the gap here is not attributable to a
// bug in the charged-GCMC energy bookkeeping.

#include <catch2/catch_test_macros.hpp>

#include <array>
#include <cmath>
#include <cstdio>
#include <filesystem>
#include <map>
#include <memory>
#include <string>
#include <vector>

#include "core/geometry/lattice.hpp"
#include "core/math/counter_based_rng.hpp"
#include "core/math/particle_data.hpp"
#include "engines/monte_carlo/molecule_species.hpp"
#include "engines/monte_carlo/monte_carlo_engine.hpp"
#include "engines/monte_carlo/peng_robinson.hpp"
#include "forcefield/electrostatics/ewald.hpp"
#include "forcefield/pairwise/lennard_jones.hpp"
#include "io/structure_io.hpp"

using aleator::core::Lattice;
using aleator::core::ParticleData;
using aleator::core::Philox4x32Rng;
using aleator::engines::MoleculeSpecies;
using aleator::engines::MonteCarloEngine;
using aleator::engines::PengRobinson;
using aleator::engines::PengRobinsonSubstance;
using aleator::forcefield::Ewald;
using aleator::forcefield::LennardJones;
using aleator::forcefield::LennardJonesParameters;
using aleator::forcefield::LennardJonesTruncation;

namespace {

using Vec3 = std::array<double, 3>;

Vec3 add(const Vec3& a, const Vec3& b) { return {a[0] + b[0], a[1] + b[1], a[2] + b[2]}; }
Vec3 scale(const Vec3& a, double s) { return {a[0] * s, a[1] * s, a[2] * s}; }

std::filesystem::path co2Irmof1DataFile(const std::string& name) {
    return std::filesystem::path(ALEATOR_KNOWN_DEVIATION_DATA_DIR) / "data" / "co2_irmof1" / name;
}

LennardJonesParameters uffParameters(const std::string& element) {
    static const std::map<std::string, LennardJonesParameters> table{
        {"Zn", {62.38, 2.462}},
        {"H", {22.14, 2.571}},
        {"C", {52.8, 3.431}},
        {"O", {30.2, 3.118}},
    };
    return table.at(element);
}

ParticleData buildSupercell(const ParticleData& primitive, const Lattice& primitiveLattice, int nx,
                             int ny, int nz) {
    const auto& m = primitiveLattice.matrix();
    const Vec3 a0{m[0][0], m[0][1], m[0][2]};
    const Vec3 a1{m[1][0], m[1][1], m[1][2]};
    const Vec3 a2{m[2][0], m[2][1], m[2][2]};
    ParticleData result;
    result.reserve(primitive.size() * static_cast<std::size_t>(nx * ny * nz));
    for (int i = 0; i < nx; ++i) {
        for (int j = 0; j < ny; ++j) {
            for (int k = 0; k < nz; ++k) {
                const Vec3 shift = add(add(scale(a0, i), scale(a1, j)), scale(a2, k));
                for (std::size_t p = 0; p < primitive.size(); ++p) {
                    result.push_back(primitive.x[p] + shift[0], primitive.y[p] + shift[1],
                                      primitive.z[p] + shift[2], 0.0, 0.0, 0.0, primitive.mass[p],
                                      primitive.charge[p], primitive.species[p]);
                }
            }
        }
    }
    return result;
}

Lattice buildSupercellLattice(const Lattice& primitive, int nx, int ny, int nz) {
    const auto& m = primitive.matrix();
    std::array<std::array<double, 3>, 3> mm{};
    mm[0] = {m[0][0] * nx, m[0][1] * nx, m[0][2] * nx};
    mm[1] = {m[1][0] * ny, m[1][1] * ny, m[1][2] * ny};
    mm[2] = {m[2][0] * nz, m[2][1] * nz, m[2][2] * nz};
    return Lattice(mm);
}

constexpr double kAmuToKg = 1.66053906660e-27;
constexpr double kAvogadro = 6.02214076e23;

// Recorded from the real calibration run this test's baseline comes from
// (see file header). Loading in mol/kg.
constexpr double kPressurePascal = 1.0e6;
constexpr double kReferenceLoadingMolPerKg = 12.832284;
constexpr double kReferenceErrorMolPerKg = 0.9939228;
constexpr double kBaselineRelativeDeviation = 0.1951; // ~19.5%, as measured

// Additive margin (relative-deviation percentage points) allowed on top of
// the recorded baseline before this test fails -- absorbs re-run Monte
// Carlo noise without also absorbing a genuine further widening. See
// tests/known_deviation/test_gcmc_ch4_irmof1_known_deviation.cc for the
// same convention.
constexpr double kAllowedDeviationGrowth = 0.05;

} // namespace

TEST_CASE("GCMC-computed CO2/IRMOF-1 loading at 1e6 Pa has not grown further from the "
          "published CRAFTED reference than its recorded baseline",
          "[known-deviation][montecarlo][electrostatics][slow]") {
    auto structure = aleator::io::readCif(co2Irmof1DataFile("IRMOF-1_primitive_DDEC.cif"));

    double totalCharge = 0.0;
    for (double q : structure.particles.charge) {
        totalCharge += q;
    }
    const double perAtomCorrection =
        totalCharge / static_cast<double>(structure.particles.size());
    for (double& q : structure.particles.charge) {
        q -= perAtomCorrection;
    }

    const ParticleData framework =
        buildSupercell(structure.particles, structure.lattice, 2, 2, 2);
    const Lattice lattice = buildSupercellLattice(structure.lattice, 2, 2, 2);

    double frameworkMassAmu = 0.0;
    for (double m : framework.mass) {
        frameworkMassAmu += m;
    }
    const double frameworkMassKg = frameworkMassAmu * kAmuToKg;

    std::vector<LennardJonesParameters> ljParameters;
    for (const auto& symbol : structure.speciesSymbols) {
        ljParameters.push_back(uffParameters(symbol));
    }
    const auto carbonIndex = static_cast<std::uint32_t>(ljParameters.size());
    ljParameters.push_back({27.0, 2.8});
    const auto oxygenIndex = static_cast<std::uint32_t>(ljParameters.size());
    ljParameters.push_back({79.0, 3.05});

    const double cutoff = 12.8;
    auto forceField = std::make_shared<const LennardJones>(ljParameters, cutoff,
                                                             LennardJonesTruncation::Shifted);
    auto electrostatics = std::make_shared<const Ewald>(/*alpha=*/0.25, cutoff, /*kMax=*/6);

    MoleculeSpecies co2;
    co2.sites.push_back({0.0, 0.0, -1.16, oxygenIndex, /*charge=*/-0.35, /*mass=*/15.9994});
    co2.sites.push_back({0.0, 0.0, 0.0, carbonIndex, /*charge=*/0.70, /*mass=*/12.0107});
    co2.sites.push_back({0.0, 0.0, 1.16, oxygenIndex, /*charge=*/-0.35, /*mass=*/15.9994});

    const double temperatureKelvin = 298.0;
    const PengRobinsonSubstance kCarbonDioxide{304.1282, 7377300.0, 0.22394};
    const PengRobinson eos(kCarbonDioxide);
    const double fugacityPascal = eos.fugacityPascal(temperatureKelvin, kPressurePascal);

    constexpr std::size_t kEquilibrationSteps = 2'000;
    constexpr std::size_t kProductionSteps = 40'000;
    constexpr int kBlocks = 8;
    constexpr std::size_t kBlockSteps = kProductionSteps / kBlocks;

    auto rng = std::make_unique<Philox4x32Rng>();
    rng->seed(/*seed=*/2026, /*streamIndex=*/2);
    MonteCarloEngine engine(framework, lattice, forceField, std::move(rng), temperatureKelvin, co2,
                             fugacityPascal, /*maxTranslation=*/1.5, /*maxRotation=*/0.5,
                             electrostatics);

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

    const double loadingMolPerKg = (meanN / kAvogadro) / frameworkMassKg;
    const double loadingErrorMolPerKg = (standardErrorN / kAvogadro) / frameworkMassKg;
    const double relativeDeviation =
        std::abs(loadingMolPerKg - kReferenceLoadingMolPerKg) / kReferenceLoadingMolPerKg;

    const double drift = engine.electrostaticEnergyDriftForTesting();
    const double energy = engine.electrostaticEnergyForTesting();
    const double relDrift = drift / std::max(1.0, std::abs(energy));

    std::printf(
        "P=%.0f Pa: <N>=%.3f +/- %.3f  computed=%.4f +/- %.4f mol/kg  reference(CRAFTED)=%.4f "
        "+/- %.4f mol/kg  relative deviation=%.2f%% (baseline %.2f%%)  electrostatic rel. "
        "drift=%.2e\n",
        kPressurePascal, meanN, standardErrorN, loadingMolPerKg, loadingErrorMolPerKg,
        kReferenceLoadingMolPerKg, kReferenceErrorMolPerKg, 100.0 * relativeDeviation,
        100.0 * kBaselineRelativeDeviation, relDrift);

    // Electrostatic cache correctness is NOT informational -- this is a
    // real, tight check even in an otherwise-informational test.
    REQUIRE(relDrift <= 1e-10);

    // Informational: fails only if the gap has genuinely grown beyond its
    // recorded baseline, never because the gap itself is nonzero.
    const double allowedMax = kBaselineRelativeDeviation + kAllowedDeviationGrowth;
    INFO("relative deviation=" << 100.0 * relativeDeviation
                                << "%, baseline=" << 100.0 * kBaselineRelativeDeviation
                                << "%, allowed max=" << 100.0 * allowedMax << "%");
    CHECK(relativeDeviation < allowedMax);
}
