#pragma once

#include <cstddef>
#include <vector>

#include "core/geometry/lattice.hpp"
#include "core/math/particle_data.hpp"
#include "forcefield/electrostatics/ewald.hpp"

namespace aleator::forcefield {

/// One charged site (position + charge), independent of core::ParticleData's
/// SoA layout or any GCMC-specific molecule type: this class is usable by
/// any single-trajectory consumer that perturbs a handful of charges at a
/// time, not just MonteCarloEngine, so it doesn't depend on engines/ types
/// (CLAUDE.md #3: core -> forcefield -> engines, never the reverse).
struct ChargedSite {
    double x = 0.0;
    double y = 0.0;
    double z = 0.0;
    double charge = 0.0;
};

/// Incrementally-maintained Ewald reciprocal-space, self-energy, and
/// intramolecular-exclusion state, across a sequence of trial moves
/// (translate/rotate one rigid molecule, insert a molecule, delete a
/// molecule) that each perturb only a handful of charged sites -- CLAUDE.md
/// hygiene-session milestone: "the reciprocal-space Ewald sum is a global
/// sum over the structure factor S(k), so it does not decompose into a
/// per-particle contribution... the standard solution is incremental
/// structure-factor bookkeeping."
///
/// SCOPE: this class owns exactly {reciprocal, self, intramolecular
/// exclusion}. It deliberately does NOT own the real-space term: that term
/// is already local/pairwise (same cost structure as LennardJones), stays
/// O(N) direct-scan via Ewald::realSpaceParticleEnergy(), and gains nothing
/// from caching. Splitting responsibility this way means each half can be
/// verified independently and neither can silently double-count the other.
///
/// WHY NOT ON Ewald ITSELF: Ewald objects are shared read-only
/// (std::shared_ptr<const ForceField>) across engine instances/replicas in
/// this codebase, and are also used as plain stateless evaluators in
/// validation tests that call computeEnergy() on many unrelated
/// configurations in sequence. A mutable per-configuration cache living on
/// Ewald itself would break both usages the moment two callers shared one
/// Ewald object. This class is instead owned by exactly one trajectory (one
/// MonteCarloEngine instance), constructed once against that trajectory's
/// starting configuration and lattice (which never change -- the framework
/// is rigid and the box is fixed for the lifetime of a GCMC run).
///
/// FRAMEWORK ASSUMPTION: the configuration passed to the constructor (the
/// framework, with zero adsorbate molecules present) contributes to
/// `committedEnergy()`'s reciprocal and self terms permanently, and
/// contributes NOTHING to the exclusion term -- MOF frameworks in a rigid
/// GCMC simulation are not modeled as "molecules" with a bonded-exclusion
/// list of their own (there is no separate bonded potential double-counting
/// anything for them); only adsorbate molecules inserted after construction
/// carry an intramolecular exclusion correction.
///
/// COMMIT/ROLLBACK SEMANTICS: every `propose*` method computes a trial
/// energy delta by calling the SAME private delta computation `commit*`
/// itself calls, and neither one ever mutates the committed cache except
/// inside `commit*`. A rejected trial therefore needs no rollback logic at
/// all: simply never call `commit*`, and the committed cache is
/// byte-for-byte whatever it was before `propose*` was called, because
/// `propose*` never touched it. This is deliberately simpler and safer than
/// mutate-then-maybe-undo -- there is no code path where an undo can be
/// forgotten or implemented slightly wrong, because there is nothing to
/// undo.
class EwaldIncrementalState {
public:
    /// Full from-scratch initialization against `particles` (typically the
    /// framework alone, zero adsorbate molecules yet) and `lattice` (fixed
    /// for the lifetime of this object). O(N * N_k).
    EwaldIncrementalState(const Ewald& ewald, const core::ParticleData& particles,
                          const core::Lattice& lattice);

    /// Reciprocal + self + exclusion energy of the current COMMITTED state.
    [[nodiscard]] double committedEnergy() const noexcept;

    /// Trial energy delta for translating/rotating one rigid molecule's
    /// sites from `oldSites` to `newSites` (same charges, same site count
    /// and order -- only positions differ). Does not mutate committed
    /// state.
    [[nodiscard]] double proposeMove(const std::vector<ChargedSite>& oldSites,
                                      const std::vector<ChargedSite>& newSites) const;

    /// Trial energy delta for inserting a new molecule at `newSites`. Does
    /// not mutate committed state -- safe to call speculatively and never
    /// commit (e.g. Widom test-particle insertion).
    [[nodiscard]] double proposeInsertion(const std::vector<ChargedSite>& newSites) const;

    /// Trial energy delta for deleting the existing molecule currently at
    /// `existingSites` (its current, committed position/charges). Does not
    /// mutate committed state.
    [[nodiscard]] double proposeDeletion(const std::vector<ChargedSite>& existingSites) const;

    /// Folds the same delta `proposeMove` would compute into the committed
    /// cache. Call only when the trial move was accepted.
    void commitMove(const std::vector<ChargedSite>& oldSites,
                     const std::vector<ChargedSite>& newSites);
    void commitInsertion(const std::vector<ChargedSite>& newSites);
    void commitDeletion(const std::vector<ChargedSite>& existingSites);

    /// Recomputes reciprocal + self + exclusion entirely from scratch
    /// against `particles` (the full current configuration) and
    /// `molecules` (current site groupings, for the exclusion term --
    /// ParticleData itself has no molecule-grouping metadata). Does NOT
    /// mutate this object. Exposed for the drift-gate validation test:
    /// compare this against committedEnergy() after many incremental
    /// updates to catch cache drift (CLAUDE.md milestone: "cache drift is
    /// exactly the kind of bug that produces smooth, plausible, wrong
    /// isotherms").
    [[nodiscard]] double recomputeFromScratch(
        const core::ParticleData& particles,
        const std::vector<std::vector<std::size_t>>& molecules) const;

    /// Replaces the committed cache with a fresh from-scratch computation
    /// against `particles`/`molecules` (i.e. recomputeFromScratch(), but
    /// mutating). Every incremental commit* is a floating-point addition;
    /// over a long enough run (this codebase's own long-run drift-gate
    /// test starts measurably exceeding 1e-10 relative error by roughly
    /// step 60,000 of a 100,000-step run, from correct-but-imperfect IEEE
    /// double accumulation, not a logic bug -- the short-run, every-step
    /// gate shows ~1e-13 agreement throughout) that accumulation needs an
    /// occasional bound. Calling this periodically (this codebase's
    /// MonteCarloEngine does so every few thousand accepted moves) is
    /// standard practice in production MD/MC codes for exactly this
    /// reason, not a workaround for a defect in the incremental math
    /// itself.
    void resync(const core::ParticleData& particles,
                const std::vector<std::vector<std::size_t>>& molecules);

    /// Raw committed S(k) (interleaved re/im) followed by self energy and
    /// exclusion energy, exposed ONLY for the rollback-exactness test
    /// (bit-for-bit comparison of the whole committed cache before and
    /// after a rejected trial) -- not intended for physics use.
    [[nodiscard]] std::vector<double> rawStateForTesting() const;

private:
    struct Delta {
        std::vector<double> deltaRe;
        std::vector<double> deltaIm;
        double reciprocalEnergy = 0.0;
        double selfEnergy = 0.0;
        double exclusionEnergy = 0.0;
        [[nodiscard]] double total() const noexcept {
            return reciprocalEnergy + selfEnergy + exclusionEnergy;
        }
    };

    /// The single source of truth for "what changes when `removed` sites
    /// (at their given, current positions) leave and `added` sites (at
    /// their given, new positions) arrive." A move is (removed=oldSites,
    /// added=newSites); an insertion is (removed={}, added=newSites); a
    /// deletion is (removed=existingSites, added={}). Every propose*/commit*
    /// pair calls this exact function, so they cannot diverge.
    [[nodiscard]] Delta computeDelta(const std::vector<ChargedSite>& removed,
                                      const std::vector<ChargedSite>& added) const;

    void apply(const Delta& delta);

    const Ewald* ewald_;
    const core::Lattice* lattice_;
    std::vector<Ewald::ReciprocalTerm> terms_;
    double volume_;
    std::vector<double> sRe_;
    std::vector<double> sIm_;
    double selfEnergy_ = 0.0;
    double exclusionEnergy_ = 0.0;
};

} // namespace aleator::forcefield
