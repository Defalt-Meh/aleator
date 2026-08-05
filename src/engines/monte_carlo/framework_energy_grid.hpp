#pragma once

#include <array>
#include <cstddef>
#include <cstdint>
#include <unordered_map>
#include <vector>

#include "core/geometry/lattice.hpp"
#include "core/math/particle_data.hpp"
#include "forcefield/pairwise/lennard_jones.hpp"

namespace aleator::engines {

/// Precomputed guest-host Lennard-Jones energy field for a rigid framework
/// (CLAUDE.md section 5 performance milestone): "the guest-host energy is a
/// fixed scalar field that can be tabulated once and interpolated... the
/// single largest available speedup" for GCMC against a host that never
/// moves. One table per guest LJ species (a multi-site adsorbate like CO2
/// needs one per site type), each a regular periodic grid of energies —
/// sum over every framework atom of pairEnergy(distance, mixedParameters)
/// — evaluated at that grid node, trilinearly interpolated at arbitrary
/// query positions in between.
///
/// Construction cost is O(gridPointCount * frameworkAtomCount) per guest
/// species — real work, meant to be paid once (at MonteCarloEngine
/// construction, or shared across several MonteCarloEngines built against
/// the same framework/force field, e.g. an isotherm's pressure points) and
/// amortized over every subsequent MC trial move, not repeated per move the
/// way the direct scan it replaces was.
///
/// Not an exact refactor (unlike core::CellList's guest-guest wiring):
/// interpolated energies differ from a direct evaluation at the same point
/// by a real, measured, documented amount that shrinks as spacingAngstrom()
/// shrinks — see the accuracy-vs-spacing measurement in
/// tests/validation/test_framework_energy_grid.cc and the default chosen in
/// engines/monte_carlo/framework_energy_grid.cc's kDefaultSpacingAngstrom
/// doc comment.
class FrameworkEnergyGrid {
public:
    /// Every stored grid node is clamped to at most this (K) before being
    /// interpolated -- see framework_energy_grid.cc's doc comment: raw LJ
    /// repulsion at a node very close to a framework atom's core can
    /// legitimately compute 1e17+ K, and interpolating that against normal
    /// neighboring nodes produces physically meaningless multi-order-of-
    /// magnitude output. At simulation temperatures (~300 K),
    /// exp(-kEnergyCapKelvin / T) is already indistinguishable from exactly
    /// zero acceptance probability, so this cannot change which trial moves
    /// get accepted -- only avoids storing/interpolating meaningless raw
    /// magnitudes. Public so callers comparing interpolated energies
    /// against a direct evaluation (e.g. this class's own accuracy test)
    /// can apply the same cap for a physically-meaningful comparison,
    /// rather than being dominated by cosmetically enormous but Boltzmann-
    /// irrelevant deviations deep in the repulsive core.
    static constexpr double kEnergyCapKelvin = 1.0e5;

    /// Builds one grid per species id in `guestSpeciesIds` (duplicates
    /// ignored). `frameworkParticles` must contain only the fixed host
    /// atoms (never adsorbate sites) -- this grid has no notion of guest
    /// occupancy, by design (a rigid framework's guest-host field doesn't
    /// depend on which/how-many guests are currently present). Throws
    /// std::invalid_argument if spacingAngstrom is not positive.
    FrameworkEnergyGrid(const forcefield::LennardJones& forceField,
                         const core::ParticleData& frameworkParticles, const core::Lattice& lattice,
                         const std::vector<std::uint32_t>& guestSpeciesIds, double spacingAngstrom);

    /// Trilinearly-interpolated guest-host energy (K) for a guest of
    /// `guestSpeciesId` at `position` (Cartesian Å, any periodic image).
    /// O(1): 8 grid-node lookups plus a weighted sum, replacing the
    /// O(frameworkAtomCount) direct scan this class exists to amortize.
    /// Throws std::out_of_range if `guestSpeciesId` wasn't included in
    /// `guestSpeciesIds` at construction.
    [[nodiscard]] double interpolate(std::uint32_t guestSpeciesId,
                                      const std::array<double, 3>& position) const;

    [[nodiscard]] double spacingAngstrom() const noexcept { return spacingAngstrom_; }
    [[nodiscard]] std::array<int, 3> pointCounts() const noexcept { return pointCounts_; }

private:
    struct Table {
        std::vector<double> values; // flattened, size nx*ny*nz, index = (iz*ny + iy)*nx + ix
    };

    [[nodiscard]] std::size_t flatten(int ix, int iy, int iz) const noexcept {
        return (static_cast<std::size_t>(iz) * static_cast<std::size_t>(pointCounts_[1]) +
                static_cast<std::size_t>(iy)) *
                   static_cast<std::size_t>(pointCounts_[0]) +
               static_cast<std::size_t>(ix);
    }

    core::Lattice lattice_;
    double spacingAngstrom_;
    std::array<int, 3> pointCounts_{1, 1, 1};
    std::unordered_map<std::uint32_t, Table> tables_;
};

} // namespace aleator::engines
