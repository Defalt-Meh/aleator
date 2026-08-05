#pragma once

#include <array>
#include <cstddef>
#include <vector>

#include "core/geometry/lattice.hpp"
#include "core/math/particle_data.hpp"

namespace aleator::core {

/// Linked-cell spatial index. Partitions the periodic cell into a grid of
/// sub-cells, sized from Lattice::perpendicularWidth (not raw lattice
/// vector length) so that any two particles within `radius` of each other
/// — under minimum image, accounting for triclinic skew — are guaranteed
/// to lie in the same or an index-adjacent sub-cell along every axis. This
/// is the same quantity Lattice::validateCutoff is built on, generalizing
/// the classic (orthorhombic) linked-cell method correctly to triclinic
/// cells (CLAUDE.md #10).
///
/// CellList only narrows candidates down; it does not compute exact
/// distances. Callers (e.g. VerletList) must still confirm each candidate
/// with an exact minimum-image distance check — forEachCandidatePair()
/// guarantees no true within-`radius` pair is ever missed (a superset
/// guarantee), not that every reported candidate is truly within `radius`.
class CellList {
public:
    /// Rebuilds the grid for `particles` under `lattice`, sized so cell
    /// extents are >= `radius` along every axis. Only particles with index
    /// >= `startIndex` are binned (default: all of them) — the returned
    /// indices are still global indices into `particles`, unshifted, so a
    /// caller can build a grid over a subrange (e.g. only the mobile/guest
    /// tail of a combined framework+adsorbate ParticleData, per CLAUDE.md
    /// section 3/5: rigid-framework GCMC's guest-guest term is the part of
    /// the energy that genuinely needs a neighbor structure — the
    /// guest-host term is handled separately, see
    /// engines/monte_carlo/monte_carlo_engine.cc) without a second
    /// ParticleData or an index-remapping step. Throws
    /// std::invalid_argument if `radius` is not positive.
    void build(const ParticleData& particles, const Lattice& lattice, double radius,
               std::size_t startIndex = 0);

    /// Invokes `callback(i, j)` — with i < j, both std::size_t particle
    /// indices — exactly once for every unordered pair of distinct
    /// particles whose sub-cells are the same or (index-)adjacent under
    /// periodic wraparound. Every pair within `radius` (build-time,
    /// minimum-image) is guaranteed to be visited; some visited pairs may
    /// be farther apart than `radius` (candidates, not confirmed
    /// neighbors).
    template <typename Callback>
    void forEachCandidatePair(Callback&& callback) const {
        const int cellCount = nx_ * ny_ * nz_;
        for (int sourceLinear = 0; sourceLinear < cellCount; ++sourceLinear) {
            const CellCoord sourceCoord = unflatten(sourceLinear);
            const auto& sourceList = bins_[static_cast<std::size_t>(sourceLinear)];
            if (sourceList.empty()) {
                continue;
            }
            for (const CellCoord& offset : offsets_) {
                const int targetLinear = flatten(wrap(sourceCoord, offset));
                const auto& targetList = bins_[static_cast<std::size_t>(targetLinear)];
                for (const std::size_t i : sourceList) {
                    for (const std::size_t j : targetList) {
                        if (j > i) {
                            callback(i, j);
                        }
                    }
                }
            }
        }
    }

    /// Invokes `callback(j)` for every binned index whose cell is the same
    /// as or (index-)adjacent to the cell containing `position` under this
    /// grid's periodic wraparound — the same candidate-superset guarantee
    /// as forEachCandidatePair(), specialized to a single arbitrary query
    /// point instead of every binned particle's own position. `position`
    /// need not correspond to any particle actually in the grid (the usual
    /// case: an MC trial move's proposed new site position). `lattice` must
    /// be the same lattice passed to the last build() call — this class
    /// does not store it, matching VerletList's build-time-lattice
    /// contract, since a caller normally already has it in hand and storing
    /// a second copy here would be redundant state to keep in sync.
    template <typename Callback>
    void forEachIndexNear(const Lattice& lattice, const std::array<double, 3>& position,
                           Callback&& callback) const {
        const CellCoord center = cellCoordFor(lattice, position);
        for (const CellCoord& offset : offsets_) {
            const int targetLinear = flatten(wrap(center, offset));
            for (const std::size_t j : bins_[static_cast<std::size_t>(targetLinear)]) {
                callback(j);
            }
        }
    }

    [[nodiscard]] int cellCountX() const noexcept { return nx_; }
    [[nodiscard]] int cellCountY() const noexcept { return ny_; }
    [[nodiscard]] int cellCountZ() const noexcept { return nz_; }

private:
    struct CellCoord {
        int x = 0;
        int y = 0;
        int z = 0;
    };

    [[nodiscard]] int flatten(CellCoord c) const noexcept { return (c.z * ny_ + c.y) * nx_ + c.x; }
    [[nodiscard]] CellCoord unflatten(int linear) const noexcept {
        const int x = linear % nx_;
        const int y = (linear / nx_) % ny_;
        const int z = linear / (nx_ * ny_);
        return {x, y, z};
    }
    [[nodiscard]] CellCoord wrap(CellCoord c, CellCoord offset) const noexcept;
    [[nodiscard]] CellCoord cellCoordFor(const Lattice& lattice,
                                          const std::array<double, 3>& position) const noexcept;

    int nx_ = 1;
    int ny_ = 1;
    int nz_ = 1;
    std::vector<std::vector<std::size_t>> bins_;
    std::vector<CellCoord> offsets_;
};

} // namespace aleator::core
