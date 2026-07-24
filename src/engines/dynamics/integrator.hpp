#pragma once

#include "core/geometry/lattice.hpp"
#include "core/math/particle_data.hpp"
#include "forcefield/force_field.hpp"

namespace aleator::engines {

/// Abstract MD integrator: advances particle positions/velocities by one
/// timestep `dt` (ps) under the given force field. Validated against NVE
/// energy drift and equipartition (CLAUDE.md #4: drift < 1e-4 k_BT/atom/ns)
/// before being trusted for anything else.
class Integrator {
public:
    virtual ~Integrator() = default;

    virtual void step(core::ParticleData& particles, const core::Lattice& lattice,
                       const forcefield::ForceField& forceField, double dtPicoseconds) = 0;
};

/// Velocity-Verlet, the standard symplectic integrator for NVE MD. Declared
/// only — throws NotImplemented until the NVE energy-drift and
/// force-vs-finite-difference validation tests exist (CLAUDE.md #4).
class VelocityVerletIntegrator final : public Integrator {
public:
    void step(core::ParticleData& particles, const core::Lattice& lattice,
              const forcefield::ForceField& forceField, double dtPicoseconds) override;
};

} // namespace aleator::engines
