#pragma once

#include <cstdint>
#include <vector>

#include "forcefield/force_field.hpp"

namespace aleator::forcefield {

/// Per-species Lennard-Jones parameters, epsilon in K, sigma in Å.
struct LennardJonesParameters {
    double epsilon = 0.0;
    double sigma = 0.0;
};

/// 12-6 Lennard-Jones pairwise force field. Gated on the NIST SRSW LJ
/// reference configurations (CLAUDE.md #4, tolerance 1e-8 relative) before
/// any of this is real — the reference test is written first and is
/// expected to fail against this stub.
class LennardJones final : public ForceField {
public:
    LennardJones(std::vector<LennardJonesParameters> perSpecies, double cutoff);

    [[nodiscard]] double computeEnergy(const core::ParticleData& particles,
                                        const core::Lattice& lattice,
                                        const core::NeighborList& neighbors) const override;

    void computeForces(const core::ParticleData& particles, const core::Lattice& lattice,
                        const core::NeighborList& neighbors, Forces& forcesOut) const override;

    [[nodiscard]] const std::vector<LennardJonesParameters>& perSpecies() const noexcept {
        return perSpecies_;
    }
    [[nodiscard]] double cutoff() const noexcept { return cutoff_; }

private:
    std::vector<LennardJonesParameters> perSpecies_;
    double cutoff_;
};

} // namespace aleator::forcefield
