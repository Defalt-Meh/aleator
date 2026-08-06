// Validation anchor for CLAUDE.md section 5's performance milestone's
// framework-energy-grid step: "Validate interpolated energies against
// direct evaluation at random points, and state the grid spacing needed to
// hold the isotherm within Monte Carlo noise. Make the spacing a config
// parameter with a documented default and a documented accuracy cost -- do
// not hide it."
//
// Ground truth is a direct sum over every framework atom (the same
// mixedParameters()/pairEnergy()/minimumImageDisplacement() calls
// FrameworkEnergyGrid's own construction uses, but evaluated live at each
// query point rather than precomputed on a grid and interpolated), computed
// independently in this file rather than by calling back into
// FrameworkEnergyGrid's own construction loop, so a shared bug in that loop
// can't hide from this test.

#include <catch2/catch_test_macros.hpp>

#include <algorithm>
#include <array>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <filesystem>
#include <map>
#include <random>
#include <stdexcept>
#include <string>
#include <vector>

#include "core/geometry/lattice.hpp"
#include "core/geometry/supercell.hpp"
#include "core/math/particle_data.hpp"
#include "engines/monte_carlo/framework_energy_grid.hpp"
#include "forcefield/pairwise/lennard_jones.hpp"
#include "io/structure_io.hpp"

using aleator::core::Lattice;
using aleator::core::ParticleData;
using aleator::core::minimumSupercellReplication;
using aleator::core::replicateSupercell;
using aleator::engines::FrameworkEnergyGrid;
using aleator::forcefield::LennardJones;
using aleator::forcefield::LennardJonesParameters;
using aleator::forcefield::LennardJonesTruncation;

namespace {

std::filesystem::path irmof1DataFile(const std::string& name) {
    return std::filesystem::path(ALEATOR_TEST_DATA_DIR) / "data" / "irmof1" / name;
}

/// Same RASPA2 GenericMOFs UFF-derived parameters as
/// tests/validation/test_gcmc_ch4_irmof1_isotherm.cc (see that file's
/// comment for provenance).
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

double directGuestHostEnergy(const LennardJones& forceField, const ParticleData& framework,
                              const Lattice& lattice, std::uint32_t guestSpecies,
                              const std::array<double, 3>& position) {
    double energy = 0.0;
    for (std::size_t f = 0; f < framework.size(); ++f) {
        const auto params = forceField.mixedParameters(guestSpecies, framework.species[f]);
        const auto d = lattice.minimumImageDisplacement(
            position, {framework.x[f], framework.y[f], framework.z[f]});
        const double r = std::sqrt(d[0] * d[0] + d[1] * d[1] + d[2] * d[2]);
        energy += forceField.pairEnergy(r, params);
    }
    // Capped the same way FrameworkEnergyGrid caps its own stored nodes
    // (see that class's kEnergyCapKelvin doc comment): without this, points
    // that land very close to a framework atom's core have a "true" direct
    // energy of 1e17+ K, and comparing that against a (necessarily capped,
    // smoothly-interpolated) grid value produces an enormous but physically
    // meaningless error dominated entirely by the repulsive-core
    // singularity rather than by anything the grid actually got wrong in a
    // way that could change GCMC acceptance.
    return std::min(energy, FrameworkEnergyGrid::kEnergyCapKelvin);
}

} // namespace

TEST_CASE("FrameworkEnergyGrid: interpolation accuracy and construction cost vs. spacing, on "
          "the real IRMOF-1/CH4 system",
          "[validation][montecarlo][slow]") {
    const auto structure = aleator::io::readCif(irmof1DataFile("IRMOF-1.cif"));
    REQUIRE(structure.particles.size() > 0);

    std::vector<LennardJonesParameters> ljParameters;
    for (const auto& symbol : structure.speciesSymbols) {
        ljParameters.push_back(uffParameters(symbol));
    }
    const auto ch4SpeciesIndex = static_cast<std::uint32_t>(ljParameters.size());
    ljParameters.push_back({158.5, 3.72}); // CH4_sp3, RASPA2 GenericMOFs

    const double cutoff = 12.0;
    const LennardJones forceField(ljParameters, cutoff, LennardJonesTruncation::Shifted);

    // Random query points, restricted to where a real GCMC trajectory would
    // actually propose insertions: the accessible (non-overlapping-with-
    // framework-core) region matters most, but insertions are proposed
    // uniformly over the WHOLE cell (see MonteCarloEngine::attemptInsertion),
    // so this test samples uniformly too rather than pre-filtering to "easy"
    // low-energy regions -- that would understate the real interpolation
    // error a GCMC run actually experiences.
    std::mt19937 queryRng(20260805);
    std::uniform_real_distribution<double> unit(0.0, 1.0);
    constexpr int kQueryPoints = 2000;
    std::vector<std::array<double, 3>> queryPoints;
    queryPoints.reserve(kQueryPoints);
    for (int i = 0; i < kQueryPoints; ++i) {
        queryPoints.push_back(structure.lattice.fractionalToCartesian(
            {unit(queryRng), unit(queryRng), unit(queryRng)}));
    }

    std::vector<double> directEnergies;
    directEnergies.reserve(queryPoints.size());
    for (const auto& p : queryPoints) {
        directEnergies.push_back(
            directGuestHostEnergy(forceField, structure.particles, structure.lattice,
                                   ch4SpeciesIndex, p));
    }

    // Candidate spacings, coarse to fine. Printed for every one so the
    // accuracy-vs-cost tradeoff that motivated the chosen default is
    // visible in the test log, not just asserted.
    const std::vector<double> candidateSpacings{1.0, 0.5, 0.3, 0.2};

    std::printf("%12s %14s %16s %16s %16s\n", "spacing(Ang)", "buildTime(ms)", "meanAbsErr(K)",
                "medianAbsErr(K)", "fracBelow1pctOfWell");
    for (double spacing : candidateSpacings) {
        const auto buildStart = std::chrono::steady_clock::now();
        const FrameworkEnergyGrid grid(forceField, structure.particles, structure.lattice,
                                        {ch4SpeciesIndex}, spacing);
        const auto buildEnd = std::chrono::steady_clock::now();
        const double buildMs =
            std::chrono::duration<double, std::milli>(buildEnd - buildStart).count();

        double sumAbsErr = 0.0;
        std::vector<double> absErrors;
        absErrors.reserve(queryPoints.size());
        for (std::size_t i = 0; i < queryPoints.size(); ++i) {
            const double interpolated = grid.interpolate(ch4SpeciesIndex, queryPoints[i]);
            const double err = std::abs(interpolated - directEnergies[i]);
            sumAbsErr += err;
            absErrors.push_back(err);
        }
        std::sort(absErrors.begin(), absErrors.end());
        const double meanAbsErr = sumAbsErr / static_cast<double>(absErrors.size());
        const double medianAbsErr = absErrors[absErrors.size() / 2];

        // "Within 1% of the LJ well depth" as an interpretable accuracy
        // yardstick: epsilon for the CH4/C (the dominant framework species
        // by count) mixed pair is sqrt(158.5*47.8562) ~= 87 K.
        const double wellDepthScaleK =
            std::sqrt(158.5 * 47.8562); // CH4 epsilon * C epsilon, geometric mean
        std::size_t belowOnePercent = 0;
        for (double e : absErrors) {
            if (e < 0.01 * wellDepthScaleK) {
                ++belowOnePercent;
            }
        }
        const double fracBelowOnePercent =
            static_cast<double>(belowOnePercent) / static_cast<double>(absErrors.size());

        std::printf("%12.2f %14.2f %16.4f %16.4f %16.4f\n", spacing, buildMs, meanAbsErr,
                    medianAbsErr, fracBelowOnePercent);
    }
}

// CLAUDE.md section 5's FrameworkEnergyGrid-serialization milestone: the single-unit-cell
// sweep above operates at cutoff ~= L_perp/2 (12.0 Ang cutoff against IRMOF-1's ~25.832 Ang
// cell), the exact pathological regime section 3 already names as one where "neither a cell
// list nor a precomputed grid have much to prune". This sweep re-runs the identical
// methodology (same force field, same cutoff, same random query protocol, same accuracy
// yardstick) on a 2x2x2 supercell instead, where the same 12.0 Ang cutoff sits at roughly
// cutoff / L_perp ~= 12.0 / 51.664 ~= 0.23 -- comfortably cutoff << L, not cutoff ~= L/2 --
// to test whether a coarser, cheaper spacing suffices there, as the milestone predicts.
TEST_CASE("FrameworkEnergyGrid: interpolation accuracy and construction cost vs. spacing, on "
          "a 2x2x2 IRMOF-1/CH4 supercell where cutoff << L",
          "[validation][montecarlo][slow]") {
    const auto unitCell = aleator::io::readCif(irmof1DataFile("IRMOF-1.cif"));
    REQUIRE(unitCell.particles.size() > 0);

    const double cutoff = 12.0;
    const auto minimumReplication = minimumSupercellReplication(unitCell.lattice, cutoff);
    // The milestone explicitly asks for 2x2x2; confirm that request actually clears the
    // minimum-image-safe replication for this cutoff and cell rather than silently using
    // whatever the minimum happens to be, so an assertion failure here means the milestone's
    // own premise ("2x2x2 supercell where cutoff << L") stopped holding for this system.
    REQUIRE(minimumReplication[0] <= 2);
    REQUIRE(minimumReplication[1] <= 2);
    REQUIRE(minimumReplication[2] <= 2);
    const auto supercell = replicateSupercell(unitCell.particles, unitCell.lattice, 2, 2, 2);

    std::vector<LennardJonesParameters> ljParameters;
    for (const auto& symbol : unitCell.speciesSymbols) {
        ljParameters.push_back(uffParameters(symbol));
    }
    const auto ch4SpeciesIndex = static_cast<std::uint32_t>(ljParameters.size());
    ljParameters.push_back({158.5, 3.72}); // CH4_sp3, RASPA2 GenericMOFs

    const LennardJones forceField(ljParameters, cutoff, LennardJonesTruncation::Shifted);

    std::mt19937 queryRng(20260805);
    std::uniform_real_distribution<double> unit(0.0, 1.0);
    constexpr int kQueryPoints = 2000;
    std::vector<std::array<double, 3>> queryPoints;
    queryPoints.reserve(kQueryPoints);
    for (int i = 0; i < kQueryPoints; ++i) {
        queryPoints.push_back(supercell.lattice.fractionalToCartesian(
            {unit(queryRng), unit(queryRng), unit(queryRng)}));
    }

    std::vector<double> directEnergies;
    directEnergies.reserve(queryPoints.size());
    for (const auto& p : queryPoints) {
        directEnergies.push_back(directGuestHostEnergy(forceField, supercell.particles,
                                                         supercell.lattice, ch4SpeciesIndex, p));
    }

    // Coarser candidates than the single-cell sweep: the milestone predicts the needed
    // spacing here is coarser and cheaper, so the candidate list is widened upward (2.0,
    // 1.5) as well as retaining the single-cell sweep's finest point (0.2) as a sanity upper
    // bound on accuracy.
    const std::vector<double> candidateSpacings{2.0, 1.5, 1.0, 0.5, 0.2};

    std::printf("%12s %14s %16s %16s %16s\n", "spacing(Ang)", "buildTime(ms)", "meanAbsErr(K)",
                "medianAbsErr(K)", "fracBelow1pctOfWell");
    for (double spacing : candidateSpacings) {
        const auto buildStart = std::chrono::steady_clock::now();
        const FrameworkEnergyGrid grid(forceField, supercell.particles, supercell.lattice,
                                        {ch4SpeciesIndex}, spacing);
        const auto buildEnd = std::chrono::steady_clock::now();
        const double buildMs =
            std::chrono::duration<double, std::milli>(buildEnd - buildStart).count();

        double sumAbsErr = 0.0;
        std::vector<double> absErrors;
        absErrors.reserve(queryPoints.size());
        for (std::size_t i = 0; i < queryPoints.size(); ++i) {
            const double interpolated = grid.interpolate(ch4SpeciesIndex, queryPoints[i]);
            const double err = std::abs(interpolated - directEnergies[i]);
            sumAbsErr += err;
            absErrors.push_back(err);
        }
        std::sort(absErrors.begin(), absErrors.end());
        const double meanAbsErr = sumAbsErr / static_cast<double>(absErrors.size());
        const double medianAbsErr = absErrors[absErrors.size() / 2];

        const double wellDepthScaleK = std::sqrt(158.5 * 47.8562);
        std::size_t belowOnePercent = 0;
        for (double e : absErrors) {
            if (e < 0.01 * wellDepthScaleK) {
                ++belowOnePercent;
            }
        }
        const double fracBelowOnePercent =
            static_cast<double>(belowOnePercent) / static_cast<double>(absErrors.size());

        std::printf("%12.2f %14.2f %16.4f %16.4f %16.4f\n", spacing, buildMs, meanAbsErr,
                    medianAbsErr, fracBelowOnePercent);
    }
}
