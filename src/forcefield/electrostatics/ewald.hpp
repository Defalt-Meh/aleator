#pragma once

#include "forcefield/force_field.hpp"

namespace aleator::forcefield {

/// Ewald summation electrostatics: real-space sum + reciprocal-space sum +
/// self-interaction correction + (if molecules have excluded intramolecular
/// pairs) exclusion correction + conducting/"tinfoil" boundary term.
///
/// CLAUDE.md #4: "Ewald is the highest-risk component in the codebase. It is
/// easy to write an implementation that is off by the self-interaction
/// term, the intramolecular exclusion correction, or the surface/tinfoil
/// boundary term — and such an implementation will produce smooth,
/// plausible, completely wrong isotherms." This class exists only as the
/// declared entry point; do not implement any of computeEnergy /
/// computeForces without first wiring the NaCl rocksalt Madelung constant
/// (1.747564594633..., 1e-6 relative tolerance) validation test, and do not
/// consider it done until that test passes — see CLAUDE.md invariant #1.
class Ewald final : public ForceField {
public:
    Ewald(double alpha, double realSpaceCutoff, int kMax);

    [[nodiscard]] double computeEnergy(const core::ParticleData& particles,
                                        const core::Lattice& lattice,
                                        const core::NeighborList& neighbors) const override;

    void computeForces(const core::ParticleData& particles, const core::Lattice& lattice,
                        const core::NeighborList& neighbors, Forces& forcesOut) const override;

    [[nodiscard]] double alpha() const noexcept { return alpha_; }
    [[nodiscard]] double realSpaceCutoff() const noexcept { return realSpaceCutoff_; }
    [[nodiscard]] int kMax() const noexcept { return kMax_; }

private:
    double alpha_;
    double realSpaceCutoff_;
    int kMax_;
};

} // namespace aleator::forcefield
