// CLAUDE.md section 5's FrameworkEnergyGrid-serialization milestone asked for the cell-list
// benchmark to be "wired to a system with real guest occupancy (CO2 at ~127 molecules), where
// it can actually show a difference": bench_gcmc_ch4_irmof1.cc's own doc comment already
// records that its guest occupancy tops out around 23 molecules, "too few for an O(N) vs.
// O(1)-ish distinction to show up against a 424-atom framework scan". This benchmark uses the
// real charged CO2/IRMOF-1 system (2x2x2 supercell, real DDEC framework charges, real Ewald
// electrostatics) at its real high-pressure point (P=1e6 Pa, see
// tests/validation/test_gcmc_co2_irmof1_isotherm.cc and CLAUDE.md's defect 1), where real
// calibration runs found ~127 molecules present at equilibrium -- the same system, same
// parameters, same seed protocol as that validation test, not an approximation.
//
// Honest limitation, not swept under the rug: unlike bench_gcmc_ch4_irmof1.cc's "before this
// milestone" profiling (a real before/after comparison), there is currently no runtime toggle
// in MonteCarloEngine to force guest-guest interactions back onto the pre-CellList O(N^2)
// scan (moleculeDispersionEnergy unconditionally uses guestCellList_ -- see that function's
// definition). So this benchmark cannot itself report a CellList-on vs. CellList-off delta;
// what it gives is a real absolute number, at real occupancy, that a future regression
// (e.g. reverting the CellList wiring) would actually move -- which the 23-molecule CH4
// system could not.

#include <array>
#include <cstdint>
#include <cstdio>
#include <filesystem>
#include <map>
#include <memory>
#include <stdexcept>
#include <string>
#include <vector>

#include <benchmark/benchmark.h>

#include "core/geometry/lattice.hpp"
#include "core/math/counter_based_rng.hpp"
#include "core/math/particle_data.hpp"
#include "engines/monte_carlo/molecule_species.hpp"
#include "engines/monte_carlo/monte_carlo_engine.hpp"
#include "engines/monte_carlo/peng_robinson.hpp"
#include "forcefield/electrostatics/ewald.hpp"
#include "forcefield/pairwise/lennard_jones.hpp"
#include "io/structure_io.hpp"

namespace {

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

using Vec3 = std::array<double, 3>;

Vec3 add(const Vec3& a, const Vec3& b) { return {a[0] + b[0], a[1] + b[1], a[2] + b[2]}; }
Vec3 scale(const Vec3& a, double s) { return {a[0] * s, a[1] * s, a[2] * s}; }

std::filesystem::path co2Irmof1DataFile(const std::string& name) {
    return std::filesystem::path(ALEATOR_CO2_IRMOF1_DATA_DIR) / name;
}

// Same CRAFTED bare-UFF parameters as test_gcmc_co2_irmof1_isotherm.cc -- see that file's
// comment for provenance (deliberately not RASPA2 GenericMOFs, which the uncharged CH4
// system uses instead).
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

void BM_GcmcRunCo2HighOccupancy(benchmark::State& state) {
    auto structure = aleator::io::readCif(co2Irmof1DataFile("IRMOF-1_primitive_DDEC.cif"));

    double totalCharge = 0.0;
    for (double q : structure.particles.charge) {
        totalCharge += q;
    }
    const double perAtomCorrection = totalCharge / static_cast<double>(structure.particles.size());
    for (double& q : structure.particles.charge) {
        q -= perAtomCorrection;
    }

    const ParticleData framework = buildSupercell(structure.particles, structure.lattice, 2, 2, 2);
    const Lattice lattice = buildSupercellLattice(structure.lattice, 2, 2, 2);

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
    // The real high-pressure point from test_gcmc_co2_irmof1_isotherm.cc's own reference
    // isotherm (CLAUDE.md defect 1): P=1e6 Pa, ~127 molecules at equilibrium -- the exact
    // occupancy regime this benchmark exists to exercise.
    const double fugacityPascal = eos.fugacityPascal(temperatureKelvin, /*pressurePascal=*/1.0e6);

    // 2000-step equilibration matches the validation test's own calibrated value (see that
    // file's comment on kEquilibrationSteps); the numSteps benchmark argument covers only
    // production, timed after equilibration completes off the clock.
    const auto numSteps = static_cast<std::size_t>(state.range(0));

    for (auto _ : state) {
        state.PauseTiming();
        auto rng = std::make_unique<Philox4x32Rng>();
        rng->seed(/*seed=*/2026, /*streamIndex=*/0);
        MonteCarloEngine engine(framework, lattice, forceField, std::move(rng), temperatureKelvin,
                                 co2, fugacityPascal, /*maxTranslation=*/1.5,
                                 /*maxRotation=*/0.5, electrostatics);
        engine.run(2'000); // equilibration, off the clock
        state.ResumeTiming();

        engine.run(numSteps);
        std::printf("[CO2/IRMOF-1 high-occupancy] molecule count after equilibration + %zu "
                    "production steps: %zu\n",
                    numSteps, engine.moleculeCount());
        benchmark::DoNotOptimize(engine.moleculeCount());
    }
    state.SetItemsProcessed(static_cast<std::int64_t>(state.iterations()) *
                             static_cast<std::int64_t>(numSteps));
}

} // namespace

// 5,000 production steps: enough to report a stable per-step cost at real ~127-molecule
// occupancy without making the benchmark suite itself slow (this system's per-step cost is
// dominated by Ewald real+reciprocal terms and a much larger guest-guest candidate count than
// the CH4 benchmark ever reaches).
BENCHMARK(BM_GcmcRunCo2HighOccupancy)
    ->Arg(5'000)
    ->Unit(benchmark::kMillisecond)
    ->Iterations(1);
