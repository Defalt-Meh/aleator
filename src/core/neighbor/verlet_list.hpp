#pragma once

#include <vector>

#include "core/geometry/lattice.hpp"
#include "core/math/particle_data.hpp"
#include "core/neighbor/neighbor_list.hpp"

namespace aleator::core {

/// Verlet (skin) neighbor list. build() uses a CellList sized to
/// cutoff + skin to find every pair within cutoff + skin (a candidate
/// superset, cached internally); pairs() filters that cached superset down
/// to the exact cutoff against the CURRENT contents of the ParticleData
/// passed to build() every time it's called (not a value cached at build
/// time) — so pairs() stays exactly correct across any number of small
/// particle moves, without needing to rebuild the superset itself, as long
/// as needsRebuild() is false.
///
/// **Aliasing contract**: build() stores a pointer to the `particles`
/// argument, not a copy. The caller must keep mutating that *same* object
/// in place (the standard pattern: one engine-owned ParticleData, moved
/// each step) until the next build() call; passing a different object to
/// needsRebuild() or destroying/reallocating the one passed to build() is
/// undefined behavior.
///
/// Rebuild trigger (CLAUDE.md #10 — displacement-based, not a fixed
/// interval): needsRebuild() returns true once any single particle has
/// moved >= skin/2 (minimum-image distance) since the last build(). This
/// is the standard, provably-safe Verlet criterion: two particles that
/// were farther apart than cutoff + skin at build time need their combined
/// approach to exceed skin before they could be closer than cutoff: if
/// every particle's own displacement is under skin/2, no pair's combined
/// approach can exceed skin, so nothing outside the build-time cutoff+skin
/// superset can have entered the cutoff sphere.
class VerletList final : public NeighborList {
public:
    /// `skinDistance` must be positive. Larger values rebuild less often
    /// at the cost of a larger (and thus more expensive to filter)
    /// candidate superset — a tuning knob for a future session; this one
    /// just needs it to be correct for any positive value.
    explicit VerletList(double skinDistance);

    void build(const ParticleData& particles, const Lattice& lattice, double cutoff) override;
    [[nodiscard]] bool needsRebuild(const ParticleData& particles) const override;
    [[nodiscard]] const std::vector<NeighborPair>& pairs() const override;

    [[nodiscard]] double skinDistance() const noexcept { return skin_; }
    [[nodiscard]] double cutoff() const noexcept { return cutoff_; }

private:
    double skin_;
    double cutoff_ = 0.0;
    bool built_ = false;
    const ParticleData* particles_ = nullptr;
    Lattice lattice_;
    std::vector<double> referenceX_;
    std::vector<double> referenceY_;
    std::vector<double> referenceZ_;
    std::vector<NeighborPair> candidatePairs_; // within cutoff + skin, as of the last build()
    mutable std::vector<NeighborPair> exactPairsCache_;
};

} // namespace aleator::core
