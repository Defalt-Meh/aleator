#pragma once

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
    /// extents are >= `radius` along every axis. Throws
    /// std::invalid_argument if `radius` is not positive.
    void build(const ParticleData& particles, const Lattice& lattice, double radius);

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

    int nx_ = 1;
    int ny_ = 1;
    int nz_ = 1;
    std::vector<std::vector<std::size_t>> bins_;
    std::vector<CellCoord> offsets_;
};

} // namespace aleator::core
