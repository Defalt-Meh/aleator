#pragma once

#include <cstddef>
#include <memory>
#include <vector>

#include "core/geometry/lattice.hpp"
#include "core/math/counter_based_rng.hpp"
#include "core/math/particle_data.hpp"
#include "engines/monte_carlo/molecule_species.hpp"
#include "forcefield/force_field.hpp"

namespace aleator::engines {

enum class Ensemble { Nvt, Npt, Gcmc };

/// Monte Carlo engine over NVT/NPT/GCMC moves. This milestone implements
/// GCMC (grand canonical: insertion, deletion, translation, rotation of a
/// single rigid adsorbate species against a fixed host framework); NVT/NPT
/// remain declared-but-unimplemented (run() throws NotImplemented for
/// those ensembles), constructed via the original (pre-GCMC) constructor
/// preserved below for that reason.
///
/// GCMC acceptance criteria (engines/monte_carlo/gcmc_acceptance.hpp) and
/// their detailed-balance property (tests/validation/
/// test_gcmc_detailed_balance.cc) are CLAUDE.md #4-gated: "MC moves
/// violating detailed balance... produces convergent, publishable-looking,
/// incorrect isotherms." The chemical-potential/fugacity term comes from
/// engines/monte_carlo/peng_robinson.hpp, not raw external pressure.
class MonteCarloEngine {
public:
    /// General (NVT/NPT placeholder) constructor — unimplemented ensembles
    /// only; run() throws NotImplemented. Preserved from the scaffolding
    /// milestone.
    MonteCarloEngine(Ensemble ensemble, core::ParticleData particles, core::Lattice lattice,
                      std::shared_ptr<const forcefield::ForceField> forceField,
                      std::unique_ptr<core::CounterBasedRng> rng, double temperatureKelvin);

    /// GCMC constructor. `frameworkParticles` are the fixed host-structure
    /// atoms — never touched by any move, always present. `species`
    /// describes the single rigid adsorbate molecule GCMC inserts/deletes
    /// (a single-site species, e.g. united-atom CH4, has exactly one site
    /// at the origin — rotation is then a formal no-op since a sphere's
    /// energy is orientation-independent, but the move is still attempted
    /// and still detailed-balance-correct). `fugacityPascal` is the
    /// pre-computed gas-phase fugacity (e.g. via PengRobinson::
    /// fugacityPascal) at the desired external pressure; converted once,
    /// here, to this codebase's internal K/Angstrom^3 energy-density
    /// convention via kPascalToInternal.
    MonteCarloEngine(core::ParticleData frameworkParticles, core::Lattice lattice,
                      std::shared_ptr<const forcefield::ForceField> forceField,
                      std::unique_ptr<core::CounterBasedRng> rng, double temperatureKelvin,
                      MoleculeSpecies species, double fugacityPascal,
                      double maxTranslationDisplacementAngstrom = 1.5,
                      double maxRotationAngleRadians = 0.5);

    /// Runs `numSteps` MC steps. For Ensemble::Gcmc, each step picks one of
    /// {translation, rotation, insertion, deletion} uniformly at random and
    /// attempts it (a no-op if there are currently zero molecules and the
    /// picked move needs one). Throws NotImplemented for Nvt/Npt.
    void run(std::size_t numSteps);

    [[nodiscard]] const core::ParticleData& particles() const noexcept { return particles_; }
    [[nodiscard]] Ensemble ensemble() const noexcept { return ensemble_; }
    [[nodiscard]] double temperatureKelvin() const noexcept { return temperatureKelvin_; }
    [[nodiscard]] std::size_t moleculeCount() const noexcept { return molecules_.size(); }
    [[nodiscard]] std::size_t frameworkParticleCount() const noexcept { return frameworkCount_; }
    [[nodiscard]] double fugacityInternal() const noexcept { return fugacityInternal_; }

    /// Pa -> this codebase's internal K/Angstrom^3 pressure/energy-density
    /// convention: divide by k_B (SI 2019 exact, 1.380649e-23 J/K) and by
    /// 1e30 (Angstrom^3 per m^3) — i.e. P[K/Ang^3] = P[Pa] / (k_B * 1e30).
    static constexpr double kPascalToInternal = 1.0 / (1.380649e-23 * 1.0e30);

    /// Zero-coverage Henry coefficient via Widom test-particle insertion
    /// against the framework ALONE (any adsorbate molecules currently
    /// present are ignored — the standard zero-coverage definition):
    ///   K_H = (1/T) * mean(exp(-deltaU/T))
    /// over `numTrials` random insertions (uniform position, and — for a
    /// multi-site species — uniform random orientation). Units
    /// molecules/(Angstrom^3 * (K/Angstrom^3)) = 1/K. Exposed for
    /// tests/validation/test_gcmc_widom_henry.cc's self-consistency check
    /// against the low-pressure isotherm slope.
    [[nodiscard]] double widomInsertionHenryCoefficient(std::size_t numTrials) const;

private:
    Ensemble ensemble_;
    core::ParticleData particles_;
    core::Lattice lattice_;
    std::shared_ptr<const forcefield::ForceField> forceField_;
    std::unique_ptr<core::CounterBasedRng> rng_;
    double temperatureKelvin_;

    // GCMC-only state (populated only by the GCMC constructor).
    std::size_t frameworkCount_ = 0;
    MoleculeSpecies species_;
    double fugacityInternal_ = 0.0;
    double maxTranslation_ = 1.0;
    double maxRotationAngle_ = 0.5;
    std::vector<std::vector<std::size_t>> molecules_;

    void attemptTranslation();
    void attemptRotation();
    void attemptInsertion();
    void attemptDeletion();
    void removeMolecule(std::size_t moleculeIndex);
};

} // namespace aleator::engines
