#pragma once

#include <cstddef>
#include <memory>

#include "core/geometry/lattice.hpp"
#include "core/math/counter_based_rng.hpp"
#include "core/math/particle_data.hpp"
#include "forcefield/force_field.hpp"

namespace aleator::engines {

enum class Ensemble { Nvt, Npt, Gcmc };

/// Monte Carlo engine over NVT/NPT/GCMC moves. Move selection, acceptance
/// criteria, and (for GCMC) insertion/deletion biasing are all
/// correctness-critical: CLAUDE.md #4 gates every move type on a
/// forward/reverse detailed-balance acceptance-ratio audit, and #10 lists
/// "MC moves violating detailed balance" as a named failure mode that
/// produces convergent, publishable-looking, *wrong* isotherms. None of
/// that exists yet — run() throws NotImplemented.
class MonteCarloEngine {
public:
    MonteCarloEngine(Ensemble ensemble, core::ParticleData particles, core::Lattice lattice,
                      std::shared_ptr<const forcefield::ForceField> forceField,
                      std::unique_ptr<core::CounterBasedRng> rng, double temperatureKelvin);

    /// Runs `numSteps` MC steps (each step attempts one move). Throws
    /// NotImplemented.
    void run(std::size_t numSteps);

    [[nodiscard]] const core::ParticleData& particles() const noexcept { return particles_; }
    [[nodiscard]] Ensemble ensemble() const noexcept { return ensemble_; }
    [[nodiscard]] double temperatureKelvin() const noexcept { return temperatureKelvin_; }

private:
    Ensemble ensemble_;
    core::ParticleData particles_;
    core::Lattice lattice_;
    std::shared_ptr<const forcefield::ForceField> forceField_;
    std::unique_ptr<core::CounterBasedRng> rng_;
    double temperatureKelvin_;
};

} // namespace aleator::engines
