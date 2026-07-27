#pragma once

#include <cstddef>
#include <vector>

#include "core/geometry/lattice.hpp"
#include "core/math/particle_data.hpp"

namespace aleator::core {

/// One (i, j) pair within cutoff, i < j. Concrete implementations decide
/// their own internal representation (cell lists, Verlet lists, tiled
/// variants — CLAUDE.md #3) but must all be able to produce this common
/// pair view for callers in forcefield/ and engines/.
struct NeighborPair {
    std::size_t i;
    std::size_t j;
};

/// Abstract neighbor-list interface. Pure virtual: there is deliberately no
/// default implementation to fall back on, so a caller can never get a
/// silently-empty or silently-stale neighbor list by accident — every
/// concrete strategy (cell list, Verlet list, ...) must be selected
/// explicitly.
///
/// Implementations own the skin-distance / rebuild-frequency decision
/// internally (CLAUDE.md #10: "Neighbor list not rebuilt often enough" is a
/// named failure mode — needs a skin-distance check, not a fixed rebuild
/// interval).
class NeighborList {
public:
    virtual ~NeighborList() = default;

    /// (Re)builds the neighbor list for the given particles and cell at the
    /// given cutoff (Å). `cutoff` must satisfy cutoff <= L/2 for every
    /// lattice direction under minimum-image convention; rejecting a larger
    /// cutoff is the caller's/implementation's responsibility, not
    /// something to discover at runtime (CLAUDE.md #10).
    virtual void build(const ParticleData& particles, const Lattice& lattice, double cutoff) = 0;

    /// Returns true if any particle has moved far enough since the last
    /// build() that a rebuild is required to stay correct within the
    /// implementation's skin distance.
    [[nodiscard]] virtual bool needsRebuild(const ParticleData& particles) const = 0;

    /// All pairs (i, j), i < j, with a minimum-image separation at or below
    /// the cutoff passed to the last build() — exact, not merely a
    /// candidate superset, as long as needsRebuild() is false. (A skin-based
    /// implementation may internally track a larger cutoff+skin candidate
    /// set to avoid rebuilding every call, but that is an implementation
    /// detail: what this method returns is always filtered to the exact
    /// cutoff. See the concrete class's docs for what "since the last
    /// build()" means in terms of which ParticleData it reads current
    /// positions from.)
    [[nodiscard]] virtual const std::vector<NeighborPair>& pairs() const = 0;
};

} // namespace aleator::core
