#include "engines/dynamics/integrator.hpp"

#include "core/exceptions.hpp"

namespace aleator::engines {

void VelocityVerletIntegrator::step(core::ParticleData&, const core::Lattice&,
                                     const forcefield::ForceField&, double) {
    throw aleator::NotImplemented("VelocityVerletIntegrator::step");
}

} // namespace aleator::engines
