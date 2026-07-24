#pragma once

#include <vector>

#include "core/geometry/lattice.hpp"
#include "core/math/particle_data.hpp"
#include "core/memory/aligned_allocator.hpp"
#include "core/neighbor/neighbor_list.hpp"

namespace aleator::forcefield {

/// Force output, SoA (CLAUDE.md invariant #2), one entry per particle, in
/// K/Å (internal energy units are K = energy/k_B, length in Å — CLAUDE.md
/// #2.4).
struct Forces {
    std::vector<double, core::AlignedAllocator<double>> fx, fy, fz;

    void resize(std::size_t n) {
        fx.resize(n);
        fy.resize(n);
        fz.resize(n);
    }
};

/// Abstract energy/force evaluator. Pure virtual, like NeighborList: every
/// concrete force field (Lennard-Jones, Ewald electrostatics, ...) must
/// implement both members explicitly rather than inheriting a default that
/// could silently return zero.
class ForceField {
public:
    virtual ~ForceField() = default;

    /// Total potential energy (K) of `particles` under this force field,
    /// given the periodic cell and a neighbor list already built at (at
    /// least) this force field's cutoff.
    [[nodiscard]] virtual double computeEnergy(const core::ParticleData& particles,
                                                const core::Lattice& lattice,
                                                const core::NeighborList& neighbors) const = 0;

    /// Per-particle forces (K/Å), accumulated into `forcesOut` (which the
    /// caller must have sized to particles.size() beforehand).
    virtual void computeForces(const core::ParticleData& particles, const core::Lattice& lattice,
                                const core::NeighborList& neighbors, Forces& forcesOut) const = 0;
};

} // namespace aleator::forcefield
