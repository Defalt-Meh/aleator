#pragma once

#include <cstddef>
#include <memory>
#include <vector>

#include "core/geometry/lattice.hpp"
#include "core/math/counter_based_rng.hpp"
#include "core/math/particle_data.hpp"
#include "core/neighbor/cell_list.hpp"
#include "engines/monte_carlo/framework_energy_grid.hpp"
#include "engines/monte_carlo/molecule_species.hpp"
#include "forcefield/electrostatics/ewald.hpp"
#include "forcefield/electrostatics/ewald_incremental_state.hpp"
#include "forcefield/force_field.hpp"

namespace aleator::engines {

enum class Ensemble { Nvt, Npt, Gcmc };

/// Monte Carlo engine over NVT/NPT/GCMC moves. This milestone implements
/// GCMC: grand canonical insertion, deletion, translation, and rotation of
/// one or more rigid adsorbate species against a fixed host framework, plus
/// identity-swap moves between species when more than one is configured
/// (CLAUDE.md's mixture-GCMC milestone — see gcmc_acceptance.hpp's
/// swapRatio() for the acceptance-ratio derivation). NVT/NPT remain
/// declared-but-unimplemented (run() throws NotImplemented for those
/// ensembles), constructed via the original (pre-GCMC) constructor
/// preserved below for that reason.
///
/// GCMC acceptance criteria (engines/monte_carlo/gcmc_acceptance.hpp) and
/// their detailed-balance property (tests/validation/
/// test_gcmc_detailed_balance.cc) are CLAUDE.md #4-gated: "MC moves
/// violating detailed balance... produces convergent, publishable-looking,
/// incorrect isotherms." The chemical-potential/fugacity term comes from
/// engines/monte_carlo/peng_robinson.hpp, not raw external pressure.
///
/// SINGLE-SPECIES BACKWARD COMPATIBILITY, load-bearing not incidental: the
/// single-species GCMC constructor below still produces bit-identical RNG
/// draw sequences (and therefore bit-identical trajectories/isotherm
/// numbers) to the pre-mixture implementation, verified by the existing
/// CH4/IRMOF-1 known-deviation test's exact `<N>` regression baseline. This
/// is achieved by every move's species-selection step being skipped
/// entirely (no RNG draw at all, not merely a deterministic draw) whenever
/// speciesList_.size() == 1 — see attemptInsertion()/attemptDeletion()'s
/// implementations. A 2+-species engine uses a different, genuinely new
/// draw sequence (species selection is real there), which is fine since no
/// prior baseline exists to preserve for that case.
class MonteCarloEngine {
public:
    /// One rigid adsorbate molecule currently present: which
    /// MoleculeSpecies template it's an instance of (index into the
    /// species list passed to the GCMC constructor) plus its site indices
    /// into `particles_`. Replaces the pre-mixture `std::vector<
    /// std::vector<std::size_t>>` (site-lists with no species tag) — see
    /// CLAUDE.md section 0's mixture-GCMC milestone.
    struct Molecule {
        std::uint32_t speciesIndex = 0;
        std::vector<std::size_t> sites;
    };

    /// General (NVT/NPT placeholder) constructor — unimplemented ensembles
    /// only; run() throws NotImplemented. Preserved from the scaffolding
    /// milestone.
    MonteCarloEngine(Ensemble ensemble, core::ParticleData particles, core::Lattice lattice,
                      std::shared_ptr<const forcefield::ForceField> forceField,
                      std::unique_ptr<core::CounterBasedRng> rng, double temperatureKelvin);

    /// Single-species GCMC constructor — a thin convenience wrapper around
    /// the mixture constructor below (delegates to it with one-element
    /// `{species}`/`{fugacityPascal}` vectors), preserved for source
    /// compatibility and because it is what makes the single-species
    /// bit-identical-RNG-sequence guarantee (see class doc comment)
    /// meaningful to state at all: a caller that only ever wants one
    /// species doesn't need to think about the vector form.
    MonteCarloEngine(core::ParticleData frameworkParticles, core::Lattice lattice,
                      std::shared_ptr<const forcefield::ForceField> forceField,
                      std::unique_ptr<core::CounterBasedRng> rng, double temperatureKelvin,
                      MoleculeSpecies species, double fugacityPascal,
                      double maxTranslationDisplacementAngstrom = 1.5,
                      double maxRotationAngleRadians = 0.5,
                      std::shared_ptr<const forcefield::Ewald> electrostatics = nullptr,
                      std::shared_ptr<const FrameworkEnergyGrid> frameworkEnergyGrid = nullptr);

    /// Mixture GCMC constructor. `frameworkParticles` are the fixed
    /// host-structure atoms — never touched by any move, always present.
    /// `speciesList` describes every rigid adsorbate species GCMC
    /// inserts/deletes/swaps (a single-site species, e.g. united-atom CH4,
    /// has exactly one site at the origin — rotation is then a formal
    /// no-op since a sphere's energy is orientation-independent, but the
    /// move is still attempted and still detailed-balance-correct).
    /// `fugacitiesPascal[i]` is species `i`'s pre-computed gas-phase
    /// fugacity (e.g. via PengRobinson::fugacityPascal for a single
    /// component, or PengRobinsonMixture::fugacitiesPascal for a real
    /// mixture) at the desired external partial pressure; converted once,
    /// here, to this codebase's internal K/Angstrom^3 energy-density
    /// convention via kPascalToInternal. `speciesList.size()` must equal
    /// `fugacitiesPascal.size()` and both must be non-empty.
    ///
    /// With exactly one species, moves are {translation, rotation,
    /// insertion, deletion} as before. With two or more, identity-swap
    /// moves (attemptSwap(): convert one existing molecule of a randomly
    /// chosen species into a randomly chosen different species, in place,
    /// with a fresh random orientation for the new species — see
    /// gcmc_acceptance.hpp's swapRatio()) are added to the move mix, which
    /// is what makes mixture composition (and therefore selectivity)
    /// equilibrate at a practical rate; pure insertion/deletion alone
    /// equilibrates composition only very slowly at high loading.
    ///
    /// `electrostatics`, if given, adds Coulomb interactions on top of
    /// `forceField`'s dispersion energy: real-space (direct O(N) scan, via
    /// Ewald::realSpaceParticleEnergy — the same cost structure as
    /// `forceField`'s own dispersion term) plus reciprocal-space/self/
    /// exclusion (via an incrementally-maintained EwaldIncrementalState —
    /// see that class's doc comment for why the reciprocal term needs a
    /// dedicated cache rather than going through `forceField`-style
    /// per-particle calls). `forceField` itself is unaffected and still
    /// independently required to support single-particle energies
    /// (invariant #10) for the dispersion term; `electrostatics` is
    /// exempt from that check because it is never routed through
    /// `ForceField::computeParticleEnergy` at all — see CLAUDE.md section
    /// 0. Every adsorbate site's charge, plus `frameworkParticles`'
    /// charges, must sum to (numerically) zero -- checked at construction,
    /// PER SPECIES individually (a mixture inserts/deletes/swaps whole
    /// molecules of one species at a time, so each species must itself be
    /// net-neutral, same reasoning as the single-species case, just looped)
    /// -- not discovered mid-run.
    ///
    /// `frameworkEnergyGrid`, if given, replaces the direct O(frameworkCount)
    /// guest-host dispersion scan with an O(1) trilinear-interpolation
    /// lookup (CLAUDE.md section 5 performance milestone) — see
    /// framework_energy_grid.hpp. It must have been built against this same
    /// `frameworkParticles`/`lattice`/`forceField` and cover every LJ
    /// species every entry in `speciesList` uses (FrameworkEnergyGrid
    /// already builds one table per guest species id, so this was
    /// mixture-ready before this milestone); unlike `electrostatics`, this
    /// is an approximation, not an exact refactor (see that class's doc
    /// comment for the measured accuracy-vs-spacing tradeoff). Building one
    /// is real, amortizable work — callers running several
    /// MonteCarloEngines against the same framework (e.g. an isotherm's
    /// pressure points) should build it once and share it, not rebuild it
    /// per engine.
    MonteCarloEngine(core::ParticleData frameworkParticles, core::Lattice lattice,
                      std::shared_ptr<const forcefield::ForceField> forceField,
                      std::unique_ptr<core::CounterBasedRng> rng, double temperatureKelvin,
                      std::vector<MoleculeSpecies> speciesList,
                      std::vector<double> fugacitiesPascal,
                      double maxTranslationDisplacementAngstrom = 1.5,
                      double maxRotationAngleRadians = 0.5,
                      std::shared_ptr<const forcefield::Ewald> electrostatics = nullptr,
                      std::shared_ptr<const FrameworkEnergyGrid> frameworkEnergyGrid = nullptr);

    /// Runs `numSteps` MC steps. For Ensemble::Gcmc: with exactly one
    /// species, each step picks one of {translation, rotation, insertion,
    /// deletion} uniformly at random (25% each, unchanged from the
    /// pre-mixture implementation — see class doc comment on bit-identical
    /// backward compatibility). With two or more species, each step picks
    /// one of {translation, rotation, insertion, deletion, identity-swap}
    /// uniformly at random (20% each). Either way, a pick is a no-op if the
    /// chosen move currently has nothing valid to act on (e.g. deletion
    /// when the chosen species has zero molecules present). Throws
    /// NotImplemented for Nvt/Npt.
    void run(std::size_t numSteps);

    [[nodiscard]] const core::ParticleData& particles() const noexcept { return particles_; }
    [[nodiscard]] Ensemble ensemble() const noexcept { return ensemble_; }
    [[nodiscard]] double temperatureKelvin() const noexcept { return temperatureKelvin_; }
    [[nodiscard]] std::size_t moleculeCount() const noexcept { return molecules_.size(); }
    [[nodiscard]] std::size_t frameworkParticleCount() const noexcept { return frameworkCount_; }

    [[nodiscard]] std::size_t speciesCount() const noexcept { return speciesList_.size(); }

    /// Number of molecules of `speciesIndex` currently present. A linear
    /// scan over `molecules_` (CLAUDE.md section 5: "never optimize
    /// speculatively" — guest counts here are O(10-150), not a measured
    /// bottleneck; an incremental per-species counter would be a second
    /// source of truth that could silently drift from `molecules_`, which
    /// is a worse failure mode than the scan's cost).
    [[nodiscard]] std::size_t moleculeCountOfSpecies(std::size_t speciesIndex) const;

    /// Pa -> this codebase's internal K/Angstrom^3 pressure/energy-density
    /// convention: divide by k_B (SI 2019 exact, 1.380649e-23 J/K) and by
    /// 1e30 (Angstrom^3 per m^3) — i.e. P[K/Ang^3] = P[Pa] / (k_B * 1e30).
    static constexpr double kPascalToInternal = 1.0 / (1.380649e-23 * 1.0e30);

    /// Zero-coverage Henry coefficient for `speciesIndex` (default 0, the
    /// only species a single-species engine has) via Widom test-particle
    /// insertion against the framework ALONE (any adsorbate molecules
    /// currently present, of any species, are ignored — the standard
    /// zero-coverage definition):
    ///   K_H = (1/T) * mean(exp(-deltaU/T))
    /// over `numTrials` random insertions (uniform position, and — for a
    /// multi-site species — uniform random orientation). Units
    /// molecules/(Angstrom^3 * (K/Angstrom^3)) = 1/K. Exposed for
    /// tests/validation/test_gcmc_widom_henry.cc's self-consistency check
    /// against the low-pressure isotherm slope, and for IAST's pure-
    /// component-isotherm fitting (CLAUDE.md mixture milestone item 5).
    [[nodiscard]] double widomInsertionHenryCoefficient(std::size_t numTrials,
                                                         std::size_t speciesIndex = 0) const;

    /// |committed electrostatic cache energy - a fresh from-scratch
    /// recomputation|, i.e. exactly the quantity the drift-gate validation
    /// test checks stays below 1e-10 relative -- 0.0 if this engine has no
    /// electrostatics. A diagnostic hook for tests, not physics API (see
    /// EwaldIncrementalState::rawStateForTesting() for the same naming
    /// convention).
    [[nodiscard]] double electrostaticEnergyDriftForTesting() const;

    /// The committed electrostatic (reciprocal + self + exclusion) cache
    /// energy itself -- 0.0 if this engine has no electrostatics. Paired
    /// with electrostaticEnergyDriftForTesting() so a test can compute a
    /// *relative* drift, matching the drift-gate's 1e-10 relative
    /// requirement rather than an arbitrary absolute one.
    [[nodiscard]] double electrostaticEnergyForTesting() const;

private:
    Ensemble ensemble_;
    core::ParticleData particles_;
    core::Lattice lattice_;
    std::shared_ptr<const forcefield::ForceField> forceField_;
    std::unique_ptr<core::CounterBasedRng> rng_;
    double temperatureKelvin_;

    // GCMC-only state (populated only by the GCMC constructor).
    std::size_t frameworkCount_ = 0;
    std::vector<MoleculeSpecies> speciesList_;
    std::vector<double> fugacitiesInternal_;
    double maxTranslation_ = 1.0;
    double maxRotationAngle_ = 0.5;
    std::vector<Molecule> molecules_;

    // Guest-guest neighbor acceleration (CLAUDE.md section 5 performance
    // milestone): forceField_'s dispersion cutoff, fetched once at
    // construction (see monte_carlo_engine.cc for why this needs a concrete
    // LennardJones rather than living on the ForceField interface).
    // frameworkIndices_ is the fixed candidate list {0, ..., frameworkCount_
    // - 1}, precomputed once since it never changes. guestCellList_ is
    // rebuilt from the current particles_/molecules_ state after every
    // accepted move (particles_ is otherwise immutable between accepts), and
    // queried per trial via CellList::forEachIndexNear -- see
    // moleculeDispersionEnergy()'s doc comment for why the candidates it
    // gathers are sorted back into ascending index order before use.
    double guestCellListRadius_ = 0.0;
    std::vector<std::size_t> frameworkIndices_;
    core::CellList guestCellList_;

    // Optional O(1) replacement for the framework-direct-scan portion of
    // moleculeDispersionEnergy() -- null means "still scan frameworkIndices_
    // directly" (backward compatible; e.g. the charged-GCMC path doesn't use
    // this yet, see CLAUDE.md section 0).
    std::shared_ptr<const FrameworkEnergyGrid> frameworkEnergyGrid_;

    // Charged-GCMC state (populated only when the GCMC constructor is
    // given a non-null `electrostatics`; ewaldState_ is null otherwise, and
    // every electrostatics-related branch is skipped).
    std::shared_ptr<const forcefield::Ewald> electrostatics_;
    std::unique_ptr<forcefield::EwaldIncrementalState> ewaldState_;
    std::size_t movesSinceEwaldResync_ = 0;

    void attemptTranslation();
    void attemptRotation();
    void attemptInsertion();
    void attemptDeletion();
    /// Identity-swap: convert one existing molecule of a randomly chosen
    /// species into a randomly chosen different species, in place. Only
    /// ever invoked by run() when speciesList_.size() >= 2. See
    /// gcmc_acceptance.hpp's swapRatio() for the acceptance-ratio
    /// derivation and the .cc for the index-bookkeeping argument (why
    /// molecules_[moleculeIndex] is reassigned to the new species BEFORE
    /// the old sites are physically removed from particles_, not after).
    void attemptSwap();
    void removeMolecule(std::size_t moleculeIndex);

    /// Physically removes `indices` (raw particle indices whose owning
    /// molecule has already been removed or reassigned) from particles_ via
    /// the same swap-with-last mechanic removeMolecule() uses, fixing up
    /// every remaining molecule's site-index references —
    /// generalized to an arbitrary raw index list not tied to one
    /// molecule's current molecules_ entry, since attemptSwap() needs to
    /// remove the OLD species' sites after molecules_[moleculeIndex] has
    /// already been reassigned to the NEW species (see attemptSwap()'s
    /// implementation comment for why that order matters). removeMolecule()
    /// is implemented in terms of this.
    void removeParticleIndices(std::vector<std::size_t> indices);

    /// All indices j such that molecules_[j].speciesIndex == speciesIndex,
    /// in ascending order. Used by attemptDeletion()/attemptSwap() to pick
    /// a molecule of a specific species uniformly at random among that
    /// species' own population (a linear scan — same "not a measured
    /// bottleneck" reasoning as moleculeCountOfSpecies()).
    [[nodiscard]] std::vector<std::size_t> moleculeIndicesOfSpecies(std::size_t speciesIndex) const;

    /// Dispersion (forceField_) energy of `sites` against everything else,
    /// via a direct scan against the framework and a CellList-accelerated
    /// scan against other guest molecules -- the hot-path replacement for
    /// the old full-O(N)-scan moleculeInteractionEnergy() (still used
    /// as-is by widomInsertionHenryCoefficient(), which has no "other
    /// guests" to accelerate against). See the .cc for the bit-identical
    /// argument. `alsoExclude`, if given, is unioned with `sites` itself in
    /// the excluded-indices set passed to the force field -- needed by
    /// attemptSwap() to evaluate the NEW species' trial energy while the
    /// OLD (about-to-vanish) molecule's sites are still physically present
    /// in particles_ but must not contribute (they won't exist in the
    /// final, post-swap state). Every other call site uses the default
    /// (self-exclusion only), unchanged from the pre-mixture behavior.
    [[nodiscard]] double moleculeDispersionEnergy(
        const std::vector<std::size_t>& sites,
        const std::vector<std::size_t>& alsoExclude = {}) const;

    /// Rebuilds guestCellList_ from the current particles_/molecules_ state.
    /// particles_ only changes when a move is accepted, so this only needs
    /// to run once per acceptance (called at the end of every accept branch
    /// in attemptTranslation/Rotation/Insertion/Deletion/Swap), plus once
    /// at construction for the initial (empty-guest) state.
    void rebuildGuestCellList();

    /// Real-space + reciprocal/self/exclusion Coulomb energy of `sites`
    /// (a rigid molecule's site indices in `particles_`) against
    /// everything else, i.e. the electrostatic analogue of
    /// moleculeInteractionEnergy() -- 0.0 if electrostatics_ is null.
    /// `alsoExclude`: see moleculeDispersionEnergy()'s doc comment above --
    /// the identical rationale applies here.
    [[nodiscard]] double moleculeElectrostaticRealSpaceEnergy(
        const std::vector<std::size_t>& sites,
        const std::vector<std::size_t>& alsoExclude = {}) const;

    /// Periodically replaces ewaldState_'s incrementally-maintained cache
    /// with a fresh from-scratch computation, bounding long-run
    /// floating-point accumulation error (EwaldIncrementalState::resync()'s
    /// doc comment). Called after every accepted move once electrostatics_
    /// is in use; a no-op when electrostatics_ is null.
    void maybeResyncEwald();
};

} // namespace aleator::engines
