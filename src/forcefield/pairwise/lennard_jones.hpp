#pragma once

#include <cstdint>
#include <vector>

#include "forcefield/force_field.hpp"

namespace aleator::forcefield {

/// Per-species Lennard-Jones parameters, epsilon in K, sigma in Å.
struct LennardJonesParameters {
    double epsilon = 0.0;
    double sigma = 0.0;

    [[nodiscard]] bool operator==(const LennardJonesParameters&) const = default;
};

/// How the potential is modified at/beyond the cutoff radius `rc`.
/// Reference behavior and naming follow the NIST SRSW Lennard-Jones
/// reference calculations (see tests/validation/test_lennard_jones_nist.cc):
enum class LennardJonesTruncation {
    /// V(r) for r <= rc, 0 beyond. Energy AND force both discontinuous at
    /// rc. Paired with an analytic long-range correction (see
    /// tailEnergyCorrection/tailVirialCorrection) to approximate the
    /// neglected tail — NIST's "LRC" scheme.
    Truncated,
    /// V(r) - V(rc) for r <= rc, 0 beyond. Energy is continuous at rc;
    /// force is identical to Truncated (subtracting a constant doesn't
    /// change the derivative) — still discontinuous at rc.
    Shifted,
    /// V(r) - V(rc) - (r - rc) V'(rc) for r <= rc, 0 beyond — the linear
    /// term is chosen so both energy AND force go continuously to zero at
    /// rc. NIST's "LFS" (linear-force-shift) scheme; conventionally no
    /// long-range correction is added on top of this (NIST reports
    /// U_LRC = 0 for it), since the shift already crudely absorbs the
    /// near-cutoff tail rather than leaving a clean, integrable remainder.
    LinearForceShifted,
};

/// How unlike-species (epsilon_i, sigma_i) and (epsilon_j, sigma_j) are
/// combined into epsilon_ij, sigma_ij.
enum class MixingRule {
    /// sigma_ij = (sigma_i + sigma_j) / 2 (Lorentz, arithmetic mean);
    /// epsilon_ij = sqrt(epsilon_i * epsilon_j) (Berthelot, geometric mean).
    /// The most common convention (e.g. OPLS, AMBER, CHARMM defaults).
    LorentzBerthelot,
    /// Both sigma_ij and epsilon_ij as the geometric mean. Used by some
    /// force fields (e.g. historical OPLS united-atom sigma combining).
    Geometric,
};

/// 12-6 Lennard-Jones pairwise force field:
/// V_LJ(r) = 4*epsilon*[(sigma/r)^12 - (sigma/r)^6].
///
/// Validated against the NIST SRSW Lennard-Jones fluid reference
/// calculations (CLAUDE.md #4, 1e-8 relative on pair energy and virial) —
/// see tests/validation/test_lennard_jones_nist.cc, which embeds the actual
/// published NIST sample configurations and reference values (fetched from
/// nist.gov; not invented). Forces are validated against a central finite
/// difference of the energy (tests/unit/test_lennard_jones_forces.cc,
/// 1e-6 tolerance).
class LennardJones final : public ForceField {
public:
    LennardJones(std::vector<LennardJonesParameters> perSpecies, double cutoff,
                 LennardJonesTruncation truncation = LennardJonesTruncation::Truncated,
                 MixingRule mixingRule = MixingRule::LorentzBerthelot);

    [[nodiscard]] std::string name() const override { return "LennardJones"; }

    /// True: computeParticleEnergy() below is a real direct O(N) pairwise
    /// sum, not the NotImplemented default.
    [[nodiscard]] bool supportsSingleParticleEnergy() const override { return true; }

    /// Sum of pairwise V(r) over every pair in `neighbors.pairs()` with
    /// minimum-image separation <= cutoff(). Does NOT include the
    /// long-range tail correction — see tailEnergyCorrection(), which
    /// callers add explicitly (matching how the NIST reference table
    /// itself reports U_pair and U_LRC as separate columns).
    [[nodiscard]] double computeEnergy(const core::ParticleData& particles,
                                        const core::Lattice& lattice,
                                        const core::NeighborList& neighbors) const override;

    void computeForces(const core::ParticleData& particles, const core::Lattice& lattice,
                        const core::NeighborList& neighbors, Forces& forcesOut) const override;

    /// Direct O(N) sum of pairEnergy(index, j) over every other particle j
    /// in `particles` (skipping `excludedIndices`) — see ForceField's doc
    /// comment for why this bypasses the neighbor list. Used by
    /// engines/monte_carlo for GCMC trial-move energy differences.
    [[nodiscard]] double computeParticleEnergy(
        std::size_t index, const core::ParticleData& particles, const core::Lattice& lattice,
        const std::vector<std::size_t>& excludedIndices = {}) const override;

    /// Same per-pair math as computeParticleEnergy (identical calls to
    /// pairEnergy()/mixedParameters()/Lattice::minimumImageDisplacement, so
    /// the two are bit-identical whenever `candidateIndices` happens to be
    /// {0, ..., particles.size()-1}\{index}), but summed only over
    /// `candidateIndices` instead of every particle. Exact as long as
    /// `candidateIndices` is a superset of every j actually within
    /// cutoff() of `index` (pairEnergy() already returns 0 beyond cutoff()
    /// for anything included that isn't) — e.g. the output of
    /// core::CellList::forEachIndexNear() built with radius >= cutoff().
    [[nodiscard]] double computeParticleEnergyOverCandidates(
        std::size_t index, const core::ParticleData& particles, const core::Lattice& lattice,
        const std::vector<std::size_t>& candidateIndices,
        const std::vector<std::size_t>& excludedIndices = {}) const override;

    /// Instantaneous pairwise virial W = -sum_{i<j} r_ij * dV/dr|_{r_ij}
    /// over `neighbors.pairs()` (same sign/definition NIST's reference
    /// table uses for W_pair). Related to the configurational pressure
    /// contribution by P_config = W / (3V).
    [[nodiscard]] double computeVirial(const core::ParticleData& particles,
                                        const core::Lattice& lattice,
                                        const core::NeighborList& neighbors) const;

    /// Analytic long-range (tail) correction to the total energy from
    /// interactions beyond the cutoff, assuming g(r) = 1 for r > cutoff():
    /// U_LRC = (2*pi/V) * sum_a sum_b N_a*N_b*epsilon_ab*sigma_ab^3 *
    ///         [(4/9)(sigma_ab/rc)^9 - (4/3)(sigma_ab/rc)^3]
    /// (sum over all ordered pairs of species a, b; reduces to the
    /// standard single-species formula (8/3)*pi*N*rho*epsilon*sigma^3*
    /// [(1/3)(sigma/rc)^9 - (sigma/rc)^3] when there is one species).
    [[nodiscard]] double tailEnergyCorrection(const core::ParticleData& particles,
                                               const core::Lattice& lattice) const;

    /// Analytic long-range correction to the virial, derived the same way
    /// as tailEnergyCorrection() but integrating r^3 * (-dV/dr) instead of
    /// r^2 * V(r):
    /// W_LRC = (2*pi/V) * sum_a sum_b N_a*N_b*epsilon_ab*sigma_ab^3 *
    ///         [(16/3)(sigma_ab/rc)^9 - 8*(sigma_ab/rc)^3]
    [[nodiscard]] double tailVirialCorrection(const core::ParticleData& particles,
                                               const core::Lattice& lattice) const;

    /// Combines per-species parameters for species `speciesI`/`speciesJ`
    /// under this object's mixing rule. Public for testing the mixing
    /// rules directly and for callers that want mixed parameters without
    /// going through a full pair loop.
    [[nodiscard]] LennardJonesParameters mixedParameters(std::uint32_t speciesI,
                                                          std::uint32_t speciesJ) const;

    /// V(r) under this object's truncation scheme, for r <= cutoff(); 0 for
    /// r > cutoff(). Exposed for direct (single-pair) testing, e.g. finite
    /// differencing against pairForceOverR().
    [[nodiscard]] double pairEnergy(double r, const LennardJonesParameters& params) const;

    /// (-dV/dr) / r under this object's truncation scheme, for r <=
    /// cutoff(); 0 for r > cutoff(). Multiplying by the (i -> j) minimum
    /// image displacement vector gives the force on i from j directly,
    /// without a separate normalization step.
    [[nodiscard]] double pairForceOverR(double r, const LennardJonesParameters& params) const;

    [[nodiscard]] const std::vector<LennardJonesParameters>& perSpecies() const noexcept {
        return perSpecies_;
    }
    [[nodiscard]] double cutoff() const noexcept { return cutoff_; }
    [[nodiscard]] LennardJonesTruncation truncation() const noexcept { return truncation_; }
    [[nodiscard]] MixingRule mixingRule() const noexcept { return mixingRule_; }

private:
    std::vector<LennardJonesParameters> perSpecies_;
    double cutoff_;
    LennardJonesTruncation truncation_;
    MixingRule mixingRule_;
};

} // namespace aleator::forcefield
