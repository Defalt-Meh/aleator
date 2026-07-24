#include "engines/monte_carlo/monte_carlo_engine.hpp"

#include <utility>

#include "core/exceptions.hpp"

namespace aleator::engines {

MonteCarloEngine::MonteCarloEngine(Ensemble ensemble, core::ParticleData particles,
                                    core::Lattice lattice,
                                    std::shared_ptr<const forcefield::ForceField> forceField,
                                    std::unique_ptr<core::CounterBasedRng> rng,
                                    double temperatureKelvin)
    : ensemble_(ensemble),
      particles_(std::move(particles)),
      lattice_(lattice),
      forceField_(std::move(forceField)),
      rng_(std::move(rng)),
      temperatureKelvin_(temperatureKelvin) {}

void MonteCarloEngine::run(std::size_t) { throw aleator::NotImplemented("MonteCarloEngine::run"); }

} // namespace aleator::engines
