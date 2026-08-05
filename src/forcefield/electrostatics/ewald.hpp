#pragma once

#include <array>
#include <cstddef>
#include <utility>
#include <vector>

#include "forcefield/force_field.hpp"

namespace aleator::forcefield {

/// Per-term decomposition of the Ewald energy, exposed so validation tests
/// can compare each physical term (e.g. against NIST SRSW's published
/// per-term SPC/E water breakdown) individually rather than only the sum —
/// far more diagnostic than a single total when something is wrong, since a
/// coefficient bug in one term won't hide behind cancellation with another.
struct EwaldEnergyBreakdown {
    double real = 0.0;       // U_real
    double reciprocal = 0.0; // U_reciprocal
    double self = 0.0;       // U_self
    double exclusion = 0.0;  // U_exclusion
    double surface = 0.0;    // U_surface (0 under tinfoil boundary conditions)
    [[nodiscard]] double total() const { return real + reciprocal + self + exclusion + surface; }
};

/// Ewald summation electrostatics for point charges under periodic
/// boundary conditions:
///
///     U = U_real + U_reciprocal + U_self + U_exclusion + U_surface
///
/// Every term is implemented explicitly (see ewald.cc — each has its own
/// clearly-labeled function) and validated first against the NaCl rocksalt
/// Madelung constant (1.747564594633..., 1e-6 relative —
/// tests/validation/test_ewald_madelung.cc) before anything else, per
/// CLAUDE.md invariant #1. CLAUDE.md #4: "Ewald is the highest-risk
/// component in the codebase... easy to write an implementation that is
/// off by the self-interaction term, the intramolecular exclusion
/// correction, or the surface/tinfoil boundary term — and such an
/// implementation will produce smooth, plausible, completely wrong
/// isotherms."
///
/// Boundary condition: tinfoil (conducting surroundings, relative
/// permittivity of the surrounding medium -> infinity), the standard
/// default for condensed-phase simulation. Under tinfoil boundary
/// conditions U_surface is exactly zero (the polarization charge on a
/// conductor at infinity exactly cancels the system's net dipole term);
/// under vacuum boundary conditions (not implemented here) it would be
/// (2*pi / (3V)) * |sum_i q_i r_i|^2. This is a deliberate, explicit
/// choice, not an omission — see computeEnergy()'s implementation.
class Ewald final : public ForceField {
public:
    /// Coulomb's constant expressed for this codebase's internal units
    /// (Å length, K = energy/k_B energy, e charge): the energy (in K)
    /// between two unit charges (in e) 1 Å apart is exactly this value.
    ///     kCoulomb = e^2 / (4*pi*eps0 * 1 Angstrom * k_B)
    /// computed from CODATA constants: e = 1.602176634e-19 C (exact, SI
    /// 2019), k_B = 1.380649e-23 J/K (exact, SI 2019), eps0 =
    /// 8.8541878128e-12 F/m (CODATA 2018; relative uncertainty ~1.5e-10,
    /// i.e. exact for any purpose this codebase has). Equivalently,
    /// e^2/(4*pi*eps0) = 14.399645... eV*Angstrom, a standard tabulated
    /// constant, converted to Kelvin via 1 eV = 11604.518... K.
    static constexpr double kCoulombConstant = 167100.9468982874;

    /// `kMax` bounds the reciprocal-lattice index vectors (n0, n1, n2)
    /// summed over: n0^2 + n1^2 + n2^2 < kMax^2 + 2, a spherical
    /// (isotropic) truncation in index space, not a per-axis cube — see
    /// ewald.cc for why (including the empirical justification for the
    /// strict inequality), and note this is what makes kMax comparable
    /// across differently-shaped cells.
    ///
    /// `exclusions` lists particle index pairs (i, j) — order doesn't
    /// matter — whose direct Coulomb interaction is excluded (e.g.
    /// intramolecular pairs already accounted for by a bonded potential).
    /// Each excluded pair must not also appear as a genuine non-bonded
    /// pair; behavior is undefined if the same (i, j) is listed twice.
    Ewald(double alpha, double realSpaceCutoff, int kMax,
          std::vector<std::pair<std::size_t, std::size_t>> exclusions = {});

    [[nodiscard]] std::string name() const override { return "Ewald"; }

    /// Does NOT override computeParticleEnergy(): the reciprocal-space term
    /// is a global sum over the structure factor of every particle and
    /// isn't decomposable into a per-particle contribution without
    /// incremental structure-factor bookkeeping this class doesn't
    /// implement (see CLAUDE.md section 0, defect 2). Leaving
    /// supportsSingleParticleEnergy() at ForceField's `false` default is
    /// the deliberate, correct choice here, not an oversight — it's what
    /// lets MonteCarloEngine's constructor reject an Ewald force field
    /// loudly instead of the engine discovering the NotImplemented hole
    /// mid-run.

    /// Throws std::invalid_argument if the system's net charge is not
    /// (numerically) zero — the reciprocal-space k=0 term is only
    /// well-defined (finite) for a neutral system; a uniform neutralizing
    /// background is not implemented.
    [[nodiscard]] double computeEnergy(const core::ParticleData& particles,
                                        const core::Lattice& lattice,
                                        const core::NeighborList& neighbors) const override;

    /// Same computation as computeEnergy(), with every term exposed
    /// separately. computeEnergy() is implemented in terms of this
    /// (returns the .total()) so the two can never disagree.
    [[nodiscard]] EwaldEnergyBreakdown computeEnergyBreakdown(
        const core::ParticleData& particles, const core::Lattice& lattice,
        const core::NeighborList& neighbors) const;

    void computeForces(const core::ParticleData& particles, const core::Lattice& lattice,
                        const core::NeighborList& neighbors, Forces& forcesOut) const override;

    [[nodiscard]] double alpha() const noexcept { return alpha_; }
    [[nodiscard]] double realSpaceCutoff() const noexcept { return realSpaceCutoff_; }
    [[nodiscard]] int kMax() const noexcept { return kMax_; }
    [[nodiscard]] const std::vector<std::pair<std::size_t, std::size_t>>& exclusions()
        const noexcept {
        return exclusions_;
    }

    /// One reciprocal-lattice vector's contribution: k itself, and
    /// exp(-k^2/(4*alpha^2))/k^2 (the weight multiplying |S(k)|^2 in
    /// U_reciprocal). Exposed publicly (rather than staying an ewald.cc
    /// implementation detail) so a persistent, incremental structure-factor
    /// cache -- necessarily owned OUTSIDE this class, see
    /// EwaldIncrementalState's doc comment for why -- can enumerate the
    /// exact same k-vectors/weights this class's own computeEnergy() uses,
    /// with no risk of the two drifting apart (e.g. if kMax's spherical-vs-
    /// cubic truncation convention is ever revisited, there is exactly one
    /// place that encodes it).
    struct ReciprocalTerm {
        std::array<double, 3> k;
        double weight;
    };

    /// The full set of reciprocal-lattice terms this object's alpha/kMax
    /// would sum over for `lattice` -- i.e. exactly what computeEnergy()
    /// builds internally, exposed for callers (see above) that need to
    /// enumerate them directly rather than going through computeEnergy().
    [[nodiscard]] std::vector<ReciprocalTerm> reciprocalTerms(
        const core::Lattice& lattice) const;

    /// Real-space (erfc(alpha*r)/r) contribution of particle `index`
    /// against every other particle in `particles` NOT in
    /// `excludedIndices`, direct O(N) scan under minimum image + this
    /// object's realSpaceCutoff() -- same contract/pattern as
    /// LennardJones::computeParticleEnergy, restricted to the real-space
    /// term only (no self energy, no exclusion correction: those are
    /// O(1)-per-molecule quantities for a rigid molecule, computed
    /// separately by EwaldIncrementalState, not per-particle here).
    ///
    /// Deliberately NOT named computeParticleEnergy() and NOT what
    /// ForceField::supportsSingleParticleEnergy() reports on: that
    /// interface's contract is "the complete energy of this particle
    /// against everything else," which Ewald's reciprocal-space term
    /// structurally cannot supply this way (CLAUDE.md section 0) --
    /// calling this alone would silently omit reciprocal, self, and
    /// exclusion contributions. It exists as an explicit, honestly-scoped
    /// helper for the one caller that already knows exactly what it is and
    /// combines it correctly with the other terms: MonteCarloEngine's
    /// charged-GCMC path.
    [[nodiscard]] double realSpaceParticleEnergy(
        std::size_t index, const core::ParticleData& particles, const core::Lattice& lattice,
        const std::vector<std::size_t>& excludedIndices = {}) const;

private:
    double alpha_;
    double realSpaceCutoff_;
    int kMax_;
    std::vector<std::pair<std::size_t, std::size_t>> exclusions_;
};

} // namespace aleator::forcefield
