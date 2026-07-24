#include "forcefield/electrostatics/ewald.hpp"

#include "core/exceptions.hpp"

namespace aleator::forcefield {

Ewald::Ewald(double alpha, double realSpaceCutoff, int kMax)
    : alpha_(alpha), realSpaceCutoff_(realSpaceCutoff), kMax_(kMax) {}

double Ewald::computeEnergy(const core::ParticleData&, const core::Lattice&,
                             const core::NeighborList&) const {
    throw aleator::NotImplemented("Ewald::computeEnergy");
}

void Ewald::computeForces(const core::ParticleData&, const core::Lattice&,
                           const core::NeighborList&, Forces&) const {
    throw aleator::NotImplemented("Ewald::computeForces");
}

} // namespace aleator::forcefield
