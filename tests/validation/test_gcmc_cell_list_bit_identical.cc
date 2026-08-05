// Validation anchor for CLAUDE.md section 5's performance milestone ("Wire
// the cell list into guest-guest energy evaluation. Assert bit-identical
// energies against the current O(N) path on the existing validation systems
// before and after -- this must be an exact refactor, not an
// approximation.").
//
// MonteCarloEngine::moleculeDispersionEnergy() (monte_carlo_engine.cc)
// replaces a full O(N) scan (LennardJones::computeParticleEnergy, which
// still exists and is still used by widomInsertionHenryCoefficient) with a
// scan restricted to the framework (direct, unchanged) plus a
// core::CellList-accelerated scan over other guest molecules. This test
// checks the exact claim directly at the LennardJones/CellList level --
// independent of MonteCarloEngine's RNG/move-acceptance plumbing -- across
// varying guest occupancy (sparse to CO2/IRMOF-1-scale ~127 molecules),
// single- and multi-site "molecules", and both cubic and triclinic cells
// (CLAUDE.md invariant #6). REQUIRE(a == b), not a tolerance: the claim is
// that LennardJones::pairEnergy() returning exactly 0.0 beyond cutoff()
// lets candidates outside cutoff be dropped from the sum with no change to
// the accumulated bit pattern (IEEE754 x + 0.0 == x exactly), as long as
// every nonzero-contributing candidate is visited in the same relative
// order as the original ascending (j = 0, ..., N-1) scan -- which is why
// moleculeDispersionEnergy() sorts CellList candidates back into ascending
// index order before use (CellList's own bin-traversal order is not index
// order).

#include <catch2/catch_test_macros.hpp>

#include <algorithm>
#include <array>
#include <numeric>
#include <random>
#include <vector>

#include "core/geometry/lattice.hpp"
#include "core/math/particle_data.hpp"
#include "core/neighbor/cell_list.hpp"
#include "forcefield/pairwise/lennard_jones.hpp"

namespace {

using aleator::core::CellList;
using aleator::core::Lattice;
using aleator::core::ParticleData;
using aleator::forcefield::LennardJones;
using aleator::forcefield::LennardJonesParameters;
using aleator::forcefield::LennardJonesTruncation;
using aleator::forcefield::MixingRule;

Lattice makeLattice(std::array<double, 3> a0, std::array<double, 3> a1, std::array<double, 3> a2) {
    std::array<std::array<double, 3>, 3> m{};
    m[0] = a0;
    m[1] = a1;
    m[2] = a2;
    return Lattice(m);
}

/// `frameworkCount` framework atoms (species 0) at indices [0,
/// frameworkCount), followed by `guestCount` single-site guest particles
/// (species 1) at indices [frameworkCount, frameworkCount + guestCount) --
/// the same layout MonteCarloEngine keeps particles_ in.
ParticleData randomSystem(std::size_t frameworkCount, std::size_t guestCount, const Lattice& lattice,
                           unsigned seed) {
    std::mt19937 rng(seed);
    std::uniform_real_distribution<double> unit(0.0, 1.0);
    ParticleData particles;
    particles.resize(frameworkCount + guestCount);
    for (std::size_t i = 0; i < particles.size(); ++i) {
        const auto cart = lattice.fractionalToCartesian({unit(rng), unit(rng), unit(rng)});
        particles.x[i] = cart[0];
        particles.y[i] = cart[1];
        particles.z[i] = cart[2];
        particles.mass[i] = 1.0;
        particles.charge[i] = 0.0;
        particles.species[i] = i < frameworkCount ? 0u : 1u;
    }
    return particles;
}

/// Computes moleculeDispersionEnergy()'s exact algorithm (framework direct
/// scan + sorted CellList-restricted guest scan, single accumulator) for
/// `sites`, independent of MonteCarloEngine.
double cellListAcceleratedEnergy(const LennardJones& lj, const ParticleData& particles,
                                  const Lattice& lattice, const CellList& guestCellList,
                                  const std::vector<std::size_t>& frameworkIndices,
                                  const std::vector<std::size_t>& sites) {
    double energy = 0.0;
    for (std::size_t s : sites) {
        std::vector<std::size_t> guestCandidates;
        guestCellList.forEachIndexNear(lattice, {particles.x[s], particles.y[s], particles.z[s]},
                                        [&](std::size_t j) { guestCandidates.push_back(j); });
        std::sort(guestCandidates.begin(), guestCandidates.end());

        std::vector<std::size_t> candidates = frameworkIndices;
        candidates.insert(candidates.end(), guestCandidates.begin(), guestCandidates.end());

        energy += lj.computeParticleEnergyOverCandidates(s, particles, lattice, candidates, sites);
    }
    return energy;
}

/// Reference: the original full-O(N)-scan sum, unchanged.
double fullScanEnergy(const LennardJones& lj, const ParticleData& particles, const Lattice& lattice,
                       const std::vector<std::size_t>& sites) {
    double energy = 0.0;
    for (std::size_t s : sites) {
        energy += lj.computeParticleEnergy(s, particles, lattice, sites);
    }
    return energy;
}

void runBitIdenticalSweep(const Lattice& lattice, const char* label) {
    // Parameters and cutoff loosely modeled on the real IRMOF-1/CH4 system
    // (CLAUDE.md validation anchor), but the physical values don't matter
    // for this test -- only that the two code paths agree exactly.
    const LennardJones lj({LennardJonesParameters{96.5, 3.8}, LennardJonesParameters{148.0, 3.73}}, 12.0,
                           LennardJonesTruncation::Shifted, MixingRule::LorentzBerthelot);
    constexpr std::size_t kFrameworkCount = 424;

    // Guest occupancy from sparse to CO2/IRMOF-1-scale (~127 molecules,
    // CLAUDE.md section 0), including a couple of degenerate corners (1, 2
    // guests).
    const std::vector<std::size_t> guestCounts{1, 2, 5, 23, 60, 127};

    for (std::size_t trialIndex = 0; trialIndex < guestCounts.size(); ++trialIndex) {
        const std::size_t guestCount = guestCounts[trialIndex];
        ParticleData particles = randomSystem(kFrameworkCount, guestCount, lattice,
                                               static_cast<unsigned>(1000 + trialIndex));

        CellList guestCellList;
        guestCellList.build(particles, lattice, lj.cutoff(), kFrameworkCount);

        std::vector<std::size_t> frameworkIndices(kFrameworkCount);
        std::iota(frameworkIndices.begin(), frameworkIndices.end(), std::size_t{0});

        // Every guest site in turn plays "the moved single-site molecule".
        for (std::size_t site = kFrameworkCount; site < particles.size(); ++site) {
            const std::vector<std::size_t> sites{site};
            const double reference = fullScanEnergy(lj, particles, lattice, sites);
            const double accelerated =
                cellListAcceleratedEnergy(lj, particles, lattice, guestCellList, frameworkIndices, sites);
            INFO(label << ": guestCount=" << guestCount << " site=" << site);
            REQUIRE(accelerated == reference); // bit-identical, not a tolerance
        }

        // A 2-site "molecule" (e.g. modeling CO2's two sites), whenever
        // there are at least 2 guest particles: excludedIndices must drop
        // BOTH sites, on both code paths, identically.
        if (guestCount >= 2) {
            const std::vector<std::size_t> pairSites{kFrameworkCount, kFrameworkCount + 1};
            const double reference = fullScanEnergy(lj, particles, lattice, pairSites);
            const double accelerated = cellListAcceleratedEnergy(lj, particles, lattice, guestCellList,
                                                                   frameworkIndices, pairSites);
            INFO(label << ": guestCount=" << guestCount << " 2-site molecule");
            REQUIRE(accelerated == reference);
        }
    }
}

} // namespace

TEST_CASE("Cell-list-accelerated guest-guest dispersion energy is bit-identical to the full "
          "O(N) scan on a cubic cell",
          "[validation][gcmc]") {
    const Lattice lattice = makeLattice({25.832, 0.0, 0.0}, {0.0, 25.832, 0.0}, {0.0, 0.0, 25.832});
    runBitIdenticalSweep(lattice, "cubic");
}

TEST_CASE("Cell-list-accelerated guest-guest dispersion energy is bit-identical to the full "
          "O(N) scan on a triclinic cell",
          "[validation][gcmc]") {
    // A deliberately skewed triclinic cell (CLAUDE.md invariant #6), same
    // rough volume scale as the cubic case above.
    const Lattice lattice = makeLattice({24.0, 0.0, 0.0}, {6.0, 22.0, 0.0}, {3.5, 4.0, 21.0});
    runBitIdenticalSweep(lattice, "triclinic");
}
