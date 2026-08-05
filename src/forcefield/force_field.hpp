#pragma once

#include <cstddef>
#include <string>
#include <vector>

#include "core/exceptions.hpp"
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

    /// Short, human-readable class name (e.g. "LennardJones", "Ewald"),
    /// used only in diagnostics — CLAUDE.md invariant #10 requires that a
    /// consumer rejecting an unsupported force field name the offending
    /// class, not just report "unsupported" with no way to tell which
    /// object was at fault.
    [[nodiscard]] virtual std::string name() const = 0;

    /// Capability query for computeParticleEnergy() AND
    /// computeParticleEnergyOverCandidates() (the same per-particle trial
    /// energy, full-scan vs. candidate-list forms — a force field that
    /// implements one is expected to implement both): a consumer (e.g.
    /// MonteCarloEngine, which needs single-particle trial-move energies)
    /// must check this at construction and reject a force field that
    /// answers false, rather than discovering the NotImplemented default
    /// below the first time a real move is attempted. CLAUDE.md invariant
    /// #10: "A base-class default that throws NotImplemented is a runtime
    /// hole disguised as a compile-time contract... the engine must reject
    /// loudly [at construction]." Defaults to false; a force field that
    /// overrides computeParticleEnergy() must also override this to return
    /// true, or callers that check the capability will never use it.
    [[nodiscard]] virtual bool supportsSingleParticleEnergy() const { return false; }

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

    /// Energy (K) of the interaction between particle `index` and every
    /// other particle in `particles`, skipping any index listed in
    /// `excludedIndices` (e.g. the other sites of the same rigid molecule
    /// as `index`, whose fixed-geometry intramolecular energy is a
    /// constant that must never enter a Monte Carlo energy difference).
    /// Deliberately independent of any NeighborList: engines/monte_carlo
    /// needs this for single-trial-molecule insertion/deletion/translation/
    /// rotation moves, where only one particle (or a handful, for a
    /// multi-site molecule) changes at a time and a full neighbor-list
    /// rebuild every trial would be wasted work for no benefit at the
    /// particle counts GCMC targets — a direct O(N) scan against every
    /// other particle is both simpler and correct, matching CLAUDE.md's
    /// performance discipline (correctness first, optimize only where a
    /// profiler points).
    ///
    /// Not implemented by every ForceField (e.g. Ewald's reciprocal-space
    /// term is a global sum over all particles' structure factor and isn't
    /// decomposable into a per-particle contribution without incremental
    /// structure-factor bookkeeping that doesn't exist yet) — the default
    /// here throws NotImplemented rather than silently returning 0
    /// (CLAUDE.md invariant #7).
    [[nodiscard]] virtual double computeParticleEnergy(
        std::size_t index, const core::ParticleData& particles, const core::Lattice& lattice,
        const std::vector<std::size_t>& excludedIndices = {}) const {
        (void)index;
        (void)particles;
        (void)lattice;
        (void)excludedIndices;
        throw aleator::NotImplemented("ForceField::computeParticleEnergy");
    }

    /// Same quantity as computeParticleEnergy (energy of `index` against
    /// every particle in `candidateIndices`, minus `excludedIndices`), but
    /// summed only over `candidateIndices` rather than every particle in
    /// `particles`. Exact — not an approximation — as long as
    /// `candidateIndices` is a superset of every particle actually within
    /// this force field's cutoff of `index` (any candidate beyond cutoff
    /// contributes exactly 0, the same as it would in the full scan): the
    /// standard neighbor-structure-accelerated form of the same sum, e.g.
    /// via core::CellList::forEachIndexNear(). Gated by the same
    /// supportsSingleParticleEnergy() capability query as
    /// computeParticleEnergy() — a force field that implements one is
    /// expected to implement both (CLAUDE.md invariant #10); default throws
    /// NotImplemented for the same reason computeParticleEnergy()'s does.
    [[nodiscard]] virtual double computeParticleEnergyOverCandidates(
        std::size_t index, const core::ParticleData& particles, const core::Lattice& lattice,
        const std::vector<std::size_t>& candidateIndices,
        const std::vector<std::size_t>& excludedIndices = {}) const {
        (void)index;
        (void)particles;
        (void)lattice;
        (void)candidateIndices;
        (void)excludedIndices;
        throw aleator::NotImplemented("ForceField::computeParticleEnergyOverCandidates");
    }
};

} // namespace aleator::forcefield
