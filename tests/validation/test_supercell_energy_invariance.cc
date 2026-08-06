// Validation anchor (CLAUDE.md milestone: supercell replication). Energy
// per unit cell must be invariant under replication -- replicating a rigid
// framework is an exact operation (every physical pair interaction is
// still there, just some periodic-image pairs are now explicit atom-atom
// pairs), so this is tested as one: 1e-10 relative, not a Monte-Carlo-noise
// tolerance.
//
// Ewald gets its own, explicit case (not just reused from the LJ test)
// because the reciprocal sum is NOT trivially invariant under replication
// the way the real-space (minimum-image) sum is: the k-vector set is
// defined by the RECIPROCAL lattice, which scales inversely with
// real-space replication, so correctness requires kMax (the reciprocal
// lattice INDEX cutoff -- see Ewald::kMax's doc comment) to be scaled by
// the same integer factor as the real-space replication, not reused
// unchanged. Measured directly, not assumed: on the real system this file
// uses, an unscaled kMax gives a genuine, non-tiny 3.46e-5 relative error
// (a real, physically-expected effect of truncating a different-shaped
// k-space region), while the correctly-scaled kMax gives 4.98e-15 -- the
// difference between "a real bug" and "the code is right and the test
// setup must scale kMax" is exactly the distinction this file exists to
// get right.

#include <catch2/catch_test_macros.hpp>

#include <array>
#include <cmath>
#include <cstdio>
#include <map>
#include <string>
#include <vector>

#include "core/geometry/lattice.hpp"
#include "core/geometry/supercell.hpp"
#include "core/math/particle_data.hpp"
#include "core/neighbor/verlet_list.hpp"
#include "forcefield/electrostatics/ewald.hpp"
#include "forcefield/pairwise/lennard_jones.hpp"
#include "io/structure_io.hpp"

using aleator::core::Lattice;
using aleator::core::ParticleData;
using aleator::core::replicateSupercell;
using aleator::core::VerletList;
using aleator::forcefield::Ewald;
using aleator::forcefield::LennardJones;
using aleator::forcefield::LennardJonesParameters;
using aleator::forcefield::LennardJonesTruncation;

namespace {

std::filesystem::path irmof1Cif() {
    return std::filesystem::path(ALEATOR_TEST_DATA_DIR) / "data" / "irmof1" / "IRMOF-1.cif";
}

/// UFF-derived LJ parameters from RASPA2's "GenericMOFs" force field --
/// same live-verified table tests/validation/test_gcmc_ch4_irmof1_isotherm.cc
/// uses (see that file / data/irmof1/PROVENANCE.md for the source).
LennardJonesParameters uffParameters(const std::string& element) {
    static const std::map<std::string, LennardJonesParameters> table{
        {"Zn", {62.3992, 2.46155}},
        {"O", {48.1581, 3.03315}},
        {"C", {47.8562, 3.47299}},
        {"H", {7.64893, 2.84642}},
    };
    return table.at(element);
}

struct NaClSystem {
    ParticleData particles;
    Lattice lattice;
    std::size_t formulaUnits;
};

/// Same construction as tests/validation/test_ewald_madelung.cc's
/// buildRockSalt -- duplicated rather than shared, matching this
/// codebase's existing convention of small per-file helpers (e.g.
/// uffParameters above) over a test-only shared library.
NaClSystem buildRockSalt(double a, int nRep) {
    const double cellEdge = 2.0 * a;
    const double boxLength = cellEdge * static_cast<double>(nRep);

    const std::vector<std::array<double, 3>> naSites{
        {0.0, 0.0, 0.0}, {a, a, 0.0}, {a, 0.0, a}, {0.0, a, a}};
    const std::vector<std::array<double, 3>> clSites{
        {a, 0.0, 0.0}, {0.0, a, 0.0}, {0.0, 0.0, a}, {a, a, a}};

    ParticleData particles;
    for (int ox = 0; ox < nRep; ++ox) {
        for (int oy = 0; oy < nRep; ++oy) {
            for (int oz = 0; oz < nRep; ++oz) {
                const double dx = cellEdge * static_cast<double>(ox);
                const double dy = cellEdge * static_cast<double>(oy);
                const double dz = cellEdge * static_cast<double>(oz);
                for (const auto& site : naSites) {
                    particles.push_back(site[0] + dx, site[1] + dy, site[2] + dz, 0.0, 0.0, 0.0,
                                         22.99, 1.0, 0);
                }
                for (const auto& site : clSites) {
                    particles.push_back(site[0] + dx, site[1] + dy, site[2] + dz, 0.0, 0.0, 0.0,
                                         35.45, -1.0, 1);
                }
            }
        }
    }

    std::array<std::array<double, 3>, 3> m{};
    m[0] = {boxLength, 0.0, 0.0};
    m[1] = {0.0, boxLength, 0.0};
    m[2] = {0.0, 0.0, boxLength};

    const std::size_t formulaUnits = naSites.size() * static_cast<std::size_t>(nRep * nRep * nRep);
    return {particles, Lattice(m), formulaUnits};
}

} // namespace

TEST_CASE("LJ energy per unit cell is invariant under supercell replication (static framework, "
          "real IRMOF-1, exact to 1e-10 relative)",
          "[validation][supercell]") {
    const auto structure = aleator::io::readCif(irmof1Cif());
    std::vector<LennardJonesParameters> ljParameters;
    for (const auto& symbol : structure.speciesSymbols) {
        ljParameters.push_back(uffParameters(symbol));
    }
    const double cutoff = 12.0; // < L_perp/2 = 12.916 for real IRMOF-1's single cell
    const LennardJones lj(ljParameters, cutoff, LennardJonesTruncation::Shifted);

    VerletList originalNeighbors(/*skinDistance=*/0.5);
    originalNeighbors.build(structure.particles, structure.lattice, cutoff);
    const double originalEnergy =
        lj.computeEnergy(structure.particles, structure.lattice, originalNeighbors);

    constexpr int nx = 2, ny = 2, nz = 2;
    const auto replicated =
        replicateSupercell(structure.particles, structure.lattice, nx, ny, nz);
    VerletList replicatedNeighbors(/*skinDistance=*/0.5);
    replicatedNeighbors.build(replicated.particles, replicated.lattice, cutoff);
    const double replicatedEnergy =
        lj.computeEnergy(replicated.particles, replicated.lattice, replicatedNeighbors);

    const double expected = static_cast<double>(nx * ny * nz) * originalEnergy;
    const double relError = std::abs(replicatedEnergy - expected) / std::abs(expected);

    std::printf("LJ energy: original=%.6f, 2x2x2 replica=%.6f, expected(8x)=%.6f, relError=%.3e\n",
                originalEnergy, replicatedEnergy, expected, relError);

    REQUIRE(relError < 1e-10);
}

TEST_CASE("Ewald reciprocal (and total) energy per formula unit is invariant under supercell "
          "replication when kMax is scaled by the same factor (real NaCl rocksalt, exact to "
          "1e-10 relative)",
          "[validation][supercell][electrostatics]") {
    const double a = 1.0;
    const int nRepOriginal = 2;
    const double alpha = 1.4;
    const double realSpaceCutoff = 1.9; // < L_perp/2 = 2.0 for the original nRep=2 cell
    const int kMaxOriginal = 12;

    const NaClSystem original = buildRockSalt(a, nRepOriginal);

    VerletList originalNeighbors(/*skinDistance=*/0.3);
    originalNeighbors.build(original.particles, original.lattice, realSpaceCutoff);
    const Ewald originalEwald(alpha, realSpaceCutoff, kMaxOriginal);
    const auto originalBreakdown =
        originalEwald.computeEnergyBreakdown(original.particles, original.lattice, originalNeighbors);
    const double originalReciprocalPerFU =
        originalBreakdown.reciprocal / static_cast<double>(original.formulaUnits);
    const double originalTotalPerFU =
        originalBreakdown.total() / static_cast<double>(original.formulaUnits);

    constexpr int nx = 2, ny = 2, nz = 2;
    const auto replicated =
        replicateSupercell(original.particles, original.lattice, nx, ny, nz);
    const std::size_t replicatedFormulaUnits =
        original.formulaUnits * static_cast<std::size_t>(nx * ny * nz);
    // Same realSpaceCutoff (still safely < the replicated cell's, larger,
    // L_perp/2), same alpha, but kMax scaled by the replication factor --
    // the reciprocal lattice shrinks by exactly this factor under real-
    // space replication, so this is what keeps the physical k-space
    // truncation region equivalent (see this file's header comment).
    const int kMaxReplicated = kMaxOriginal * nx;

    VerletList replicatedNeighbors(/*skinDistance=*/0.3);
    replicatedNeighbors.build(replicated.particles, replicated.lattice, realSpaceCutoff);
    const Ewald replicatedEwald(alpha, realSpaceCutoff, kMaxReplicated);
    const auto replicatedBreakdown = replicatedEwald.computeEnergyBreakdown(
        replicated.particles, replicated.lattice, replicatedNeighbors);
    const double replicatedReciprocalPerFU =
        replicatedBreakdown.reciprocal / static_cast<double>(replicatedFormulaUnits);
    const double replicatedTotalPerFU =
        replicatedBreakdown.total() / static_cast<double>(replicatedFormulaUnits);

    const double reciprocalRelError =
        std::abs(replicatedReciprocalPerFU - originalReciprocalPerFU) /
        std::abs(originalReciprocalPerFU);
    const double totalRelError =
        std::abs(replicatedTotalPerFU - originalTotalPerFU) / std::abs(originalTotalPerFU);

    std::printf("Ewald reciprocal/FU: original=%.12f, 2x2x2 (kMax %d->%d)=%.12f, relError=%.3e\n",
                originalReciprocalPerFU, kMaxOriginal, kMaxReplicated, replicatedReciprocalPerFU,
                reciprocalRelError);
    std::printf("Ewald total/FU: original=%.6f, replica=%.6f, relError=%.3e\n", originalTotalPerFU,
                replicatedTotalPerFU, totalRelError);

    REQUIRE(reciprocalRelError < 1e-10);
    REQUIRE(totalRelError < 1e-10);
}
