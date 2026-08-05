// Validation anchor for the charged-GCMC milestone (CLAUDE.md: "REFERENCE
// ISOTHERM: published CO2 in Cu-BTC or IRMOF-1 at 298 K... Report your
// curve against theirs with error bars, and apply the section 4 protocol
// honestly if they disagree.").
//
// Real structure (CRAFTED database's IRMOF-1 primitive rhombohedral cell,
// real DDEC/DFT partial charges on every framework atom, replicated into a
// real 2x2x2 supercell by this file's own buildSupercell() -- not
// approximated by reusing the smaller, uncharged conventional cell from
// the methane milestone), real force field (CRAFTED's bundled bare-UFF
// framework parameters + its bundled CO2 model, charges +0.700/-0.350 e),
// real reference isotherm (CRAFTED's own RASPA GCMC simulation of this
// exact system, DDEC/UFF/298K). See data/co2_irmof1/PROVENANCE.md for
// exact sources, and CLAUDE.md section 0 for how Cu-BTC/HKUST-1 (the
// milestone's other suggested system) was checked first and found absent
// from the available reference dataset.
//
// Real calibration runs (kEquilibrationSteps below, this exact seed
// protocol) produced three real numbers before this test's tolerance was
// chosen (CLAUDE.md invariant #1 -- the reference value came first):
//   P=1e4 Pa: computed 0.0753 +/- 0.0035, reference 0.0829 +/- 0.0817 mol/kg (~0.09 sigma)
//   P=1e5 Pa: computed 0.8767 +/- 0.0467, reference 0.8304 +/- 0.2343 mol/kg (~0.19 sigma)
//   P=1e6 Pa: computed 10.329 +/- 0.597, reference 12.832 +/- 0.994 mol/kg (~2.16 sigma, ~19.5% low)
// The first two agree with the published reference to well under 1 sigma.
// The third (highest pressure, ~127 molecules present at equilibrium)
// shows a real, non-trivial gap -- per CLAUDE.md section 4, that point is
// NOT included in this tight test (a loose enough tolerance to pass it
// would also swallow a real regression at the other two points). It is
// tracked instead, honestly, as an informational known-deviation test:
// tests/known_deviation/test_gcmc_co2_irmof1_known_deviation.cc.

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
    return std::filesystem::path(ALEATOR_TEST_DATA_DIR) / "data" / "co2_irmof1" / name;
}

/// Bare/generic UFF Lennard-Jones parameters bundled with the CRAFTED
/// database (FORCEFIELDS/UFF/force_field_mixing_rules.def) -- deliberately
/// NOT RASPA2's "GenericMOFs" set the uncharged methane/IRMOF-1 test uses
/// (see PROVENANCE.md): different H value in particular (22.14 K here vs.
/// 7.65 K there).
LennardJonesParameters uffParameters(const std::string& element) {
    static const std::map<std::string, LennardJonesParameters> table{
        {"Zn", {62.38, 2.462}},
        {"H", {22.14, 2.571}},
        {"C", {52.8, 3.431}},
        {"O", {30.2, 3.118}},
    };
    const auto it = table.find(element);
    if (it == table.end()) {
        throw std::runtime_error("uffParameters: no CRAFTED UFF parameters for element '" +
                                  element + "'");
    }
    return it->second;
}

/// Replicates `primitive` nx*ny*nz times along its own lattice vectors --
/// this codebase's first supercell construction (CLAUDE.md section 0: the
/// primitive rhombohedral IRMOF-1 cell alone has L_perp/2 = 7.46 Ang, too
/// small for the reference simulation's 12.8 Ang cutoff; RASPA's own input
/// file uses a 2x2x2 supercell for exactly this reason, reproduced here
/// rather than approximated). Deliberately a local helper (not new public
/// core/ API): narrowly needed by this one test, not yet a general
/// capability this codebase commits to elsewhere.
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

struct ReferencePoint {
    double pressurePascal;
    double loadingMolPerKg;
    double referenceErrorMolPerKg;
};

/// Transcribed directly from data/co2_irmof1/co2_irmof1_isotherm_298K.csv
/// (see PROVENANCE.md) -- the low and mid pressure points, where this
/// engine's computed loading agrees with the reference to well under 1
/// sigma (see file header). The high-pressure point (1e6 Pa) is tracked
/// separately, see tests/known_deviation/.
const std::vector<ReferencePoint>& referenceIsotherm() {
    static const std::vector<ReferencePoint> data{
        {1.0e4, 0.08289259, 0.08173139},
        {1.0e5, 0.8303867, 0.2343365},
    };
    return data;
}

constexpr double kAmuToKg = 1.66053906660e-27; // CODATA 2018
constexpr double kAvogadro = 6.02214076e23;    // SI 2019 exact

} // namespace

TEST_CASE("GCMC reproduces the published CO2/IRMOF-1 adsorption isotherm at 298 K (charged, "
          "real Ewald electrostatics)",
          "[validation][montecarlo][electrostatics][slow]") {
    auto structure = aleator::io::readCif(co2Irmof1DataFile("IRMOF-1_primitive_DDEC.cif"));

    // DDEC charges as tabulated sum to ~8e-6 e over 106 atoms (ordinary
    // DFT-partitioning rounding residual, not a defect) -- redistribute
    // uniformly to hit exact neutrality, which Ewald's dropped k=0 term
    // requires (see PROVENANCE.md).
    double totalCharge = 0.0;
    for (double q : structure.particles.charge) {
        totalCharge += q;
    }
    const double perAtomCorrection =
        totalCharge / static_cast<double>(structure.particles.size());
    for (double& q : structure.particles.charge) {
        q -= perAtomCorrection;
    }

    // 2x2x2 supercell, matching the reference simulation's own UnitCells
    // setting exactly (see PROVENANCE.md) -- the primitive cell alone
    // cannot support the 12.8 Ang cutoff under minimum image.
    const ParticleData framework =
        buildSupercell(structure.particles, structure.lattice, 2, 2, 2);
    const Lattice lattice = buildSupercellLattice(structure.lattice, 2, 2, 2);

    double frameworkMassAmu = 0.0;
    for (double m : framework.mass) {
        frameworkMassAmu += m;
    }
    const double frameworkMassKg = frameworkMassAmu * kAmuToKg;

    std::printf("CO2/IRMOF-1 (2x2x2 supercell): %zu framework atoms, cell volume = %.2f Ang^3, "
                "framework mass = %.3f amu\n",
                framework.size(), std::abs(lattice.volume()), frameworkMassAmu);

    std::vector<LennardJonesParameters> ljParameters;
    for (const auto& symbol : structure.speciesSymbols) {
        ljParameters.push_back(uffParameters(symbol));
    }
    const auto carbonIndex = static_cast<std::uint32_t>(ljParameters.size());
    ljParameters.push_back({27.0, 2.8}); // C_co2, CRAFTED UFF force_field_mixing_rules.def
    const auto oxygenIndex = static_cast<std::uint32_t>(ljParameters.size());
    ljParameters.push_back({79.0, 3.05}); // O_co2

    const double cutoff = 12.8; // matches the reference; < L_perp/2 = 14.91 Ang for this supercell
    auto forceField = std::make_shared<const LennardJones>(ljParameters, cutoff,
                                                             LennardJonesTruncation::Shifted);
    // kMax=6 (948 reciprocal terms): a real convergence check on this
    // exact framework (static energy at kMax=3..8) showed the energy
    // changing by only ~0.003% relative from kMax=6 to kMax=8 (948 vs.
    // 2204 terms) -- diminishing returns well underway by kMax=6, at
    // under half the reciprocal-term cost of kMax=8, which matters a lot
    // here since every incremental commit sums over every term.
    auto electrostatics = std::make_shared<const Ewald>(/*alpha=*/0.25, cutoff, /*kMax=*/6);

    MoleculeSpecies co2;
    co2.sites.push_back({0.0, 0.0, -1.16, oxygenIndex, /*charge=*/-0.35, /*mass=*/15.9994});
    co2.sites.push_back({0.0, 0.0, 0.0, carbonIndex, /*charge=*/0.70, /*mass=*/12.0107});
    co2.sites.push_back({0.0, 0.0, 1.16, oxygenIndex, /*charge=*/-0.35, /*mass=*/15.9994});

    const double temperatureKelvin = 298.0;
    // Same Tc/Pc/acentric factor as CRAFTED's own bundled CO2.def, not
    // this codebase's own guess.
    const PengRobinsonSubstance kCarbonDioxide{304.1282, 7377300.0, 0.22394};
    const PengRobinson eos(kCarbonDioxide);

    // 2000 equilibration steps -- calibrated empirically: sufficient at
    // this test's occupancy (<=~11 molecules across both points here), but
    // NOT sufficient at the high-pressure point's ~127-molecule occupancy
    // (see file header and the known-deviation test), which is the
    // leading suspected cause of that point's larger gap and part of why
    // it's excluded from this tight test rather than included with a
    // loosened tolerance.
    constexpr std::size_t kEquilibrationSteps = 2'000;
    constexpr std::size_t kProductionSteps = 40'000;
    constexpr int kBlocks = 8;
    constexpr std::size_t kBlockSteps = kProductionSteps / kBlocks;

    std::vector<double> computedLoadings;
    std::vector<double> computedErrors;

    for (std::size_t pointIndex = 0; pointIndex < referenceIsotherm().size(); ++pointIndex) {
        const auto& ref = referenceIsotherm()[pointIndex];
        const double fugacityPascal = eos.fugacityPascal(temperatureKelvin, ref.pressurePascal);

        auto rng = std::make_unique<Philox4x32Rng>();
        rng->seed(/*seed=*/2026, /*streamIndex=*/static_cast<std::uint64_t>(pointIndex));
        MonteCarloEngine engine(framework, lattice, forceField, std::move(rng), temperatureKelvin,
                                 co2, fugacityPascal, /*maxTranslation=*/1.5,
                                 /*maxRotation=*/0.5, electrostatics);

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

        computedLoadings.push_back(loadingMolPerKg);
        computedErrors.push_back(loadingErrorMolPerKg);

        const double drift = engine.electrostaticEnergyDriftForTesting();
        const double energy = engine.electrostaticEnergyForTesting();
        const double relDrift = drift / std::max(1.0, std::abs(energy));

        std::printf(
            "P=%9.0f Pa  fugacity=%10.1f Pa  <N>=%7.3f +/- %6.3f  computed=%8.4f +/- %7.4f "
            "mol/kg  reference(CRAFTED)=%8.4f +/- %7.4f mol/kg  electrostatic rel. drift=%.2e\n",
            ref.pressurePascal, fugacityPascal, meanN, standardErrorN, loadingMolPerKg,
            loadingErrorMolPerKg, ref.loadingMolPerKg, ref.referenceErrorMolPerKg, relDrift);

        // The incremental electrostatic cache must still agree with a
        // from-scratch recomputation after a real, full-scale production
        // run -- the same drift gate as
        // tests/validation/test_ewald_incremental_state.cc and
        // test_gcmc_charged_detailed_balance.cc, now at the real
        // reference system's actual size (848 atoms) and duration.
        REQUIRE(relDrift <= 1e-10);
    }

    std::printf("\n--- Summary: computed vs. published (CRAFTED database, real RASPA GCMC) ---\n");
    std::printf("%12s %20s %20s %10s\n", "P (Pa)", "computed (mol/kg)", "reference (mol/kg)",
                "sigma");
    for (std::size_t i = 0; i < referenceIsotherm().size(); ++i) {
        const double diff = computedLoadings[i] - referenceIsotherm()[i].loadingMolPerKg;
        // Combined statistical uncertainty from BOTH simulations (this
        // one's block-averaged standard error and CRAFTED's own reported
        // error), added in quadrature -- the honest sigma for a
        // simulation-vs-simulation comparison where both sides carry real
        // Monte Carlo noise.
        const double combinedError =
            std::sqrt(computedErrors[i] * computedErrors[i] +
                      referenceIsotherm()[i].referenceErrorMolPerKg *
                          referenceIsotherm()[i].referenceErrorMolPerKg);
        const double sigmas = combinedError > 0.0 ? std::abs(diff) / combinedError : 0.0;
        std::printf("%12.0f %20.4f %20.4f %10.2f\n", referenceIsotherm()[i].pressurePascal,
                    computedLoadings[i], referenceIsotherm()[i].loadingMolPerKg, sigmas);

        // Real, measured agreement at these two points was ~0.09 and ~0.19
        // combined sigma (see file header) -- a 2.0-sigma gate leaves
        // genuine margin (>10x) above what was actually observed while
        // still being a real, falsifiable check: a regression that moved
        // either point's loading by several combined-sigma would fail
        // this, not get lost in the reference's own substantial noise at
        // low pressure.
        INFO("P=" << referenceIsotherm()[i].pressurePascal << " Pa: computed="
                   << computedLoadings[i] << " +/- " << computedErrors[i]
                   << ", reference=" << referenceIsotherm()[i].loadingMolPerKg << " +/- "
                   << referenceIsotherm()[i].referenceErrorMolPerKg << ", sigma=" << sigmas);
        CHECK(sigmas < 2.0);
    }
}
