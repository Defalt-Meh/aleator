// THE GATE (CLAUDE.md #1, #4): before anything else in forcefield/
// electrostatics is trusted, this must reproduce the NaCl rocksalt
// Madelung constant to 1e-6 relative. Per the session instructions: if
// this fails, stop and debug — do not tune the tolerance, do not proceed
// to other Ewald tests, do not report success.
//
// Structure: conventional rock-salt cubic cell (edge = 2a, a = Na-Cl
// nearest-neighbor distance), replicated 3x3x3 for comfortable real-space
// cutoff room relative to the cell (a bare 1x1x1 cell only allows a cutoff
// up to a itself, which is too tight to cleanly separate real- and
// reciprocal-space convergence). The Madelung constant is an intensive,
// per-ion-pair property of the infinite lattice, so it must come out
// identical regardless of supercell size for a correct implementation —
// this is itself part of what's being checked.
//
// Definition (standard, e.g. Kittel / Ashcroft & Mermin): for N formula
// units (N cation-anion pairs), U_total = -N * M * e^2 / (4*pi*eps0*a).
// In this codebase's internal units, U_total_K = kCoulombConstant *
// U_total_natural with U_total_natural = -N*M/a, so
// M = -(U_total_K / kCoulombConstant) * a / N.

#include <catch2/catch_test_macros.hpp>

#include <array>
#include <cmath>
#include <cstdio>
#include <vector>

#include "core/geometry/lattice.hpp"
#include "core/math/particle_data.hpp"
#include "core/neighbor/verlet_list.hpp"
#include "forcefield/electrostatics/ewald.hpp"

namespace {

using aleator::core::Lattice;
using aleator::core::ParticleData;
using aleator::core::VerletList;
using aleator::forcefield::Ewald;

constexpr double kPublishedMadelungConstant = 1.747564594633;

/// Builds an nRep x nRep x nRep supercell of the rock-salt structure with
/// Na-Cl nearest-neighbor distance `a`. Returns (particles, lattice,
/// numberOfFormulaUnits).
struct NaClSystem {
    ParticleData particles;
    Lattice lattice;
    std::size_t formulaUnits;
};

NaClSystem buildRockSalt(double a, int nRep) {
    const double cellEdge = 2.0 * a;
    const double boxLength = cellEdge * static_cast<double>(nRep);

    // FCC motif within one conventional cell: 4 Na+ at (0,0,0)-type sites,
    // 4 Cl- at (a,0,0)-type sites (each Na has 6 Cl nearest neighbors at
    // distance a once periodic images are included).
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

double computeMadelungConstant(double a, int nRep, double alpha, double realSpaceCutoff,
                                int kMax) {
    const NaClSystem system = buildRockSalt(a, nRep);

    VerletList neighbors(/*skinDistance=*/0.3);
    neighbors.build(system.particles, system.lattice, realSpaceCutoff);

    const Ewald ewald(alpha, realSpaceCutoff, kMax);
    const double totalEnergyK = ewald.computeEnergy(system.particles, system.lattice, neighbors);

    const double totalEnergyNatural = totalEnergyK / Ewald::kCoulombConstant;
    return -totalEnergyNatural * a / static_cast<double>(system.formulaUnits);
}

} // namespace

TEST_CASE("Ewald reproduces the NaCl rocksalt Madelung constant to 1e-6 relative — THE GATE",
          "[validation][electrostatics][gate]") {
    const double a = 1.0; // nearest-neighbor distance, Angstrom
    const int nRep = 3;   // 3x3x3 supercell of the conventional cell
    const double alpha = 1.4;
    const double realSpaceCutoff = 2.9; // < L_perp/2 = 3*2*a/2 = 3.0
    const int kMax = 12;

    const double computed = computeMadelungConstant(a, nRep, alpha, realSpaceCutoff, kMax);

    const double relError = std::abs(computed - kPublishedMadelungConstant) /
                             kPublishedMadelungConstant;

    // Print unconditionally (not just on failure) — this is the gate, and
    // the actual numbers must be visible regardless of pass/fail.
    std::printf("Madelung constant: computed = %.12f, published = %.12f, relative error = %.3e\n",
                computed, kPublishedMadelungConstant, relError);

    REQUIRE(relError < 1e-6);
}
