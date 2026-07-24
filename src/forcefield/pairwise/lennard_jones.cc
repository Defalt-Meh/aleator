#include "forcefield/pairwise/lennard_jones.hpp"

#include <utility>

#include "core/exceptions.hpp"

namespace aleator::forcefield {

LennardJones::LennardJones(std::vector<LennardJonesParameters> perSpecies, double cutoff)
    : perSpecies_(std::move(perSpecies)), cutoff_(cutoff) {}

double LennardJones::computeEnergy(const core::ParticleData&, const core::Lattice&,
                                    const core::NeighborList&) const {
    throw aleator::NotImplemented("LennardJones::computeEnergy");
}

void LennardJones::computeForces(const core::ParticleData&, const core::Lattice&,
                                  const core::NeighborList&, Forces&) const {
    throw aleator::NotImplemented("LennardJones::computeForces");
}

} // namespace aleator::forcefield
