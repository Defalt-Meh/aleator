#include "engines/monte_carlo/monte_carlo_engine.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <numbers>
#include <stdexcept>
#include <utility>

#include "core/exceptions.hpp"
#include "engines/monte_carlo/gcmc_acceptance.hpp"

namespace aleator::engines {

namespace {

using Vec3 = std::array<double, 3>;

Vec3 add(const Vec3& a, const Vec3& b) { return {a[0] + b[0], a[1] + b[1], a[2] + b[2]}; }
Vec3 sub(const Vec3& a, const Vec3& b) { return {a[0] - b[0], a[1] - b[1], a[2] - b[2]}; }
Vec3 scale(const Vec3& a, double s) { return {a[0] * s, a[1] * s, a[2] * s}; }
double dot(const Vec3& a, const Vec3& b) { return a[0] * b[0] + a[1] * b[1] + a[2] * b[2]; }
Vec3 cross(const Vec3& a, const Vec3& b) {
    return {a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0]};
}
Vec3 position(const core::ParticleData& particles, std::size_t i) {
    return {particles.x[i], particles.y[i], particles.z[i]};
}

void setPosition(core::ParticleData& particles, std::size_t i, const Vec3& p) {
    particles.x[i] = p[0];
    particles.y[i] = p[1];
    particles.z[i] = p[2];
}

/// Rodrigues' rotation formula: rotates v by angle theta (radians) about
/// the unit vector axis.
Vec3 rotateAboutAxis(const Vec3& v, const Vec3& axis, double theta) {
    const double c = std::cos(theta);
    const double s = std::sin(theta);
    return add(add(scale(v, c), scale(cross(axis, v), s)), scale(axis, dot(axis, v) * (1.0 - c)));
}

struct Quaternion {
    double w, x, y, z;
};

/// Uniformly-distributed unit quaternion from 3 independent Uniform(0,1)
/// draws (K. Shoemake, "Uniform Random Rotations," Graphics Gems III,
/// 1992): with r1=sqrt(1-u0), r2=sqrt(u0), any assignment of
/// {r1*sin(theta1), r1*cos(theta1), r2*sin(theta2), r2*cos(theta2)} to the
/// four quaternion components is a unit quaternion (r1^2+r2^2 = 1) whose
/// induced rotation is Haar-uniform over SO(3).
Quaternion randomUnitQuaternion(double u0, double u1, double u2) {
    const double r1 = std::sqrt(1.0 - u0);
    const double r2 = std::sqrt(u0);
    const double theta1 = 2.0 * std::numbers::pi * u1;
    const double theta2 = 2.0 * std::numbers::pi * u2;
    return {r2 * std::cos(theta2), r1 * std::sin(theta1), r1 * std::cos(theta1),
            r2 * std::sin(theta2)};
}

/// Rotates v by unit quaternion q: v' = v + 2w(qv x v) + 2(qv x (qv x v)),
/// qv = (x,y,z) — the standard matrix-free quaternion rotation identity,
/// equivalent to v' = q v q^-1 for unit q.
Vec3 rotateByQuaternion(const Vec3& v, const Quaternion& q) {
    const Vec3 qv{q.x, q.y, q.z};
    const Vec3 t = scale(cross(qv, v), 2.0);
    return add(add(v, scale(t, q.w)), cross(qv, t));
}

Vec3 moleculeCentroid(const core::ParticleData& particles, const std::vector<std::size_t>& sites) {
    Vec3 c{0.0, 0.0, 0.0};
    for (std::size_t s : sites) {
        c = add(c, position(particles, s));
    }
    return scale(c, 1.0 / static_cast<double>(sites.size()));
}

/// Energy of `sites` (a set of particle indices belonging to one rigid
/// molecule) interacting with everything else in `particles`, excluding
/// interactions among `sites` themselves (their fixed intramolecular
/// geometry never changes and must never enter an MC energy difference).
double moleculeInteractionEnergy(const forcefield::ForceField& forceField,
                                  const core::ParticleData& particles, const core::Lattice& lattice,
                                  const std::vector<std::size_t>& sites) {
    double energy = 0.0;
    for (std::size_t s : sites) {
        energy += forceField.computeParticleEnergy(s, particles, lattice, sites);
    }
    return energy;
}

/// Real-space-only electrostatic analogue of moleculeInteractionEnergy(),
/// via Ewald::realSpaceParticleEnergy (see that method's doc comment for
/// why it's deliberately not routed through ForceField::computeParticleEnergy).
double electrostaticRealSpaceInteractionEnergy(const forcefield::Ewald& ewald,
                                                const core::ParticleData& particles,
                                                const core::Lattice& lattice,
                                                const std::vector<std::size_t>& sites) {
    double energy = 0.0;
    for (std::size_t s : sites) {
        energy += ewald.realSpaceParticleEnergy(s, particles, lattice, sites);
    }
    return energy;
}

std::vector<forcefield::ChargedSite> chargedSitesOf(const core::ParticleData& particles,
                                                      const std::vector<std::size_t>& sites) {
    std::vector<forcefield::ChargedSite> chargedSites;
    chargedSites.reserve(sites.size());
    for (std::size_t s : sites) {
        chargedSites.push_back(
            {particles.x[s], particles.y[s], particles.z[s], particles.charge[s]});
    }
    return chargedSites;
}

double totalCharge(const core::ParticleData& particles) {
    double total = 0.0;
    for (double q : particles.charge) {
        total += q;
    }
    return total;
}

} // namespace

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

MonteCarloEngine::MonteCarloEngine(core::ParticleData frameworkParticles, core::Lattice lattice,
                                    std::shared_ptr<const forcefield::ForceField> forceField,
                                    std::unique_ptr<core::CounterBasedRng> rng,
                                    double temperatureKelvin, MoleculeSpecies species,
                                    double fugacityPascal, double maxTranslationDisplacementAngstrom,
                                    double maxRotationAngleRadians,
                                    std::shared_ptr<const forcefield::Ewald> electrostatics)
    : ensemble_(Ensemble::Gcmc),
      particles_(std::move(frameworkParticles)),
      lattice_(lattice),
      forceField_(std::move(forceField)),
      rng_(std::move(rng)),
      temperatureKelvin_(temperatureKelvin),
      frameworkCount_(particles_.size()),
      species_(std::move(species)),
      fugacityInternal_(fugacityPascal * kPascalToInternal),
      maxTranslation_(maxTranslationDisplacementAngstrom),
      maxRotationAngle_(maxRotationAngleRadians),
      electrostatics_(std::move(electrostatics)) {
    if (species_.sites.empty()) {
        throw std::invalid_argument("MonteCarloEngine: GCMC species must have at least one site");
    }
    if (!(temperatureKelvin_ > 0.0)) {
        throw std::invalid_argument("MonteCarloEngine: temperatureKelvin must be positive");
    }
    if (!(fugacityInternal_ > 0.0)) {
        throw std::invalid_argument("MonteCarloEngine: fugacityPascal must be positive");
    }
    // CLAUDE.md invariant #10: every trial move below (translation,
    // rotation, insertion, deletion) calls forceField_->computeParticleEnergy()
    // for the moved molecule against the framework. A force field that
    // doesn't support that (e.g. Ewald -- see ewald.hpp) must be rejected
    // here, at construction, not discovered as a NotImplemented throw the
    // first time run() is called -- "discovering an unsupported combination
    // six hours into a run is a defect."
    if (!forceField_->supportsSingleParticleEnergy()) {
        throw std::invalid_argument(
            "MonteCarloEngine: force field '" + forceField_->name() +
            "' does not support single-particle trial-move energies "
            "(supportsSingleParticleEnergy() == false), which GCMC's "
            "translation/rotation/insertion/deletion moves all require. "
            "This is not yet implemented for this force field -- see "
            "CLAUDE.md section 0.");
    }
    if (electrostatics_ != nullptr) {
        // Every adsorbate site's charge must sum to (numerically) zero:
        // GCMC inserts/deletes whole molecules one at a time, and if a
        // single molecule carried net charge, every insertion would drift
        // the system's total charge away from neutrality, which Ewald's
        // reciprocal-space k=0 term (silently dropped, see ewald.cc) is
        // only well-defined for. Caught here, once, rather than as a
        // slowly-compounding physics error deep into a run.
        double speciesCharge = 0.0;
        for (const auto& site : species_.sites) {
            speciesCharge += site.charge;
        }
        constexpr double kNeutralityTolerance = 1e-8;
        if (std::abs(speciesCharge) > kNeutralityTolerance) {
            throw std::invalid_argument(
                "MonteCarloEngine: GCMC adsorbate species is not individually charge-neutral "
                "(total site charge = " +
                std::to_string(speciesCharge) +
                " e); every insertion/deletion would drift the system away from "
                "neutrality, which electrostatics requires");
        }
        const double initialCharge = totalCharge(particles_);
        if (std::abs(initialCharge) > kNeutralityTolerance) {
            throw std::invalid_argument(
                "MonteCarloEngine: framework is not charge-neutral (total charge = " +
                std::to_string(initialCharge) + " e)");
        }
        ewaldState_ =
            std::make_unique<forcefield::EwaldIncrementalState>(*electrostatics_, particles_, lattice_);
    }
}

double MonteCarloEngine::moleculeElectrostaticRealSpaceEnergy(
    const std::vector<std::size_t>& sites) const {
    if (electrostatics_ == nullptr) {
        return 0.0;
    }
    return electrostaticRealSpaceInteractionEnergy(*electrostatics_, particles_, lattice_, sites);
}

void MonteCarloEngine::maybeResyncEwald() {
    if (ewaldState_ == nullptr) {
        return;
    }
    // Bounds long-run floating-point accumulation in the incrementally-
    // maintained reciprocal/self/exclusion cache -- see
    // EwaldIncrementalState::resync()'s doc comment. 500 was enough to
    // keep tests/validation/test_ewald_incremental_state.cc's small
    // synthetic system (~80 charges) under the 1e-10 relative gate, but
    // measurably was not enough for the real CO2/IRMOF-1 system (~850-880
    // charges depending on occupancy): a real calibration run measured
    // ~9.3e-10 relative drift with a 500-move interval at that scale.
    // More charges means more terms summed per reciprocal-term update, so
    // a larger system needs a shorter interval for the same bound -- 100
    // gives real margin under 1e-10 at this system's scale, confirmed by
    // tests/validation/test_gcmc_co2_irmof1_isotherm.cc's own drift
    // assertion.
    constexpr std::size_t kResyncInterval = 100;
    ++movesSinceEwaldResync_;
    if (movesSinceEwaldResync_ >= kResyncInterval) {
        ewaldState_->resync(particles_, molecules_);
        movesSinceEwaldResync_ = 0;
    }
}

void MonteCarloEngine::run(std::size_t numSteps) {
    if (ensemble_ != Ensemble::Gcmc) {
        throw aleator::NotImplemented("MonteCarloEngine::run (NVT/NPT)");
    }
    for (std::size_t step = 0; step < numSteps; ++step) {
        const double pick = rng_->nextUniform();
        if (pick < 0.25) {
            attemptTranslation();
        } else if (pick < 0.5) {
            attemptRotation();
        } else if (pick < 0.75) {
            attemptInsertion();
        } else {
            attemptDeletion();
        }
    }
}

void MonteCarloEngine::attemptTranslation() {
    if (molecules_.empty()) {
        return;
    }
    const auto moleculeIndex =
        static_cast<std::size_t>(rng_->nextUniform() * static_cast<double>(molecules_.size()));
    const auto& sites = molecules_[moleculeIndex];

    const double uBefore = moleculeInteractionEnergy(*forceField_, particles_, lattice_, sites);
    const double electroBefore = moleculeElectrostaticRealSpaceEnergy(sites);
    const auto oldChargedSites =
        ewaldState_ != nullptr ? chargedSitesOf(particles_, sites) : std::vector<forcefield::ChargedSite>{};

    const Vec3 delta{(2.0 * rng_->nextUniform() - 1.0) * maxTranslation_,
                      (2.0 * rng_->nextUniform() - 1.0) * maxTranslation_,
                      (2.0 * rng_->nextUniform() - 1.0) * maxTranslation_};
    std::vector<Vec3> original;
    original.reserve(sites.size());
    for (std::size_t s : sites) {
        original.push_back(position(particles_, s));
        setPosition(particles_, s, add(original.back(), delta));
    }

    const double uAfter = moleculeInteractionEnergy(*forceField_, particles_, lattice_, sites);
    const double electroAfter = moleculeElectrostaticRealSpaceEnergy(sites);
    const auto newChargedSites =
        ewaldState_ != nullptr ? chargedSitesOf(particles_, sites) : std::vector<forcefield::ChargedSite>{};
    const double reciprocalDelta =
        ewaldState_ != nullptr ? ewaldState_->proposeMove(oldChargedSites, newChargedSites) : 0.0;
    const double deltaU = (uAfter - uBefore) + (electroAfter - electroBefore) + reciprocalDelta;

    const double ratio = gcmc::metropolisRatio(temperatureKelvin_, deltaU);
    if (rng_->nextUniform() < gcmc::clampToProbability(ratio)) {
        if (ewaldState_ != nullptr) {
            ewaldState_->commitMove(oldChargedSites, newChargedSites);
            maybeResyncEwald();
        }
        return; // accepted, keep the new positions
    }
    for (std::size_t k = 0; k < sites.size(); ++k) {
        setPosition(particles_, sites[k], original[k]);
    }
}

void MonteCarloEngine::attemptRotation() {
    if (molecules_.empty() || species_.sites.size() < 2) {
        return; // a single-site species' energy is orientation-independent
    }
    const auto moleculeIndex =
        static_cast<std::size_t>(rng_->nextUniform() * static_cast<double>(molecules_.size()));
    const auto& sites = molecules_[moleculeIndex];

    const double uBefore = moleculeInteractionEnergy(*forceField_, particles_, lattice_, sites);
    const double electroBefore = moleculeElectrostaticRealSpaceEnergy(sites);
    const auto oldChargedSites =
        ewaldState_ != nullptr ? chargedSitesOf(particles_, sites) : std::vector<forcefield::ChargedSite>{};

    const Vec3 centroid = moleculeCentroid(particles_, sites);
    // Uniformly random axis (via the same quaternion-sphere-point trick
    // restricted to the vector part would bias the distribution, so pick
    // the axis directly via the standard rejection-free spherical
    // parameterization instead) and a small angle in
    // [-maxRotationAngle_, maxRotationAngle_] -- symmetric under negation,
    // which is what makes this move's proposal density symmetric and
    // Metropolis (not Metropolis-Hastings) acceptance valid.
    const double cosTheta = 2.0 * rng_->nextUniform() - 1.0;
    const double sinTheta = std::sqrt(std::max(0.0, 1.0 - cosTheta * cosTheta));
    const double phi = 2.0 * std::numbers::pi * rng_->nextUniform();
    const Vec3 axis{sinTheta * std::cos(phi), sinTheta * std::sin(phi), cosTheta};
    const double angle = (2.0 * rng_->nextUniform() - 1.0) * maxRotationAngle_;

    std::vector<Vec3> original;
    original.reserve(sites.size());
    for (std::size_t s : sites) {
        const Vec3 p = position(particles_, s);
        original.push_back(p);
        setPosition(particles_, s, add(centroid, rotateAboutAxis(sub(p, centroid), axis, angle)));
    }

    const double uAfter = moleculeInteractionEnergy(*forceField_, particles_, lattice_, sites);
    const double electroAfter = moleculeElectrostaticRealSpaceEnergy(sites);
    const auto newChargedSites =
        ewaldState_ != nullptr ? chargedSitesOf(particles_, sites) : std::vector<forcefield::ChargedSite>{};
    const double reciprocalDelta =
        ewaldState_ != nullptr ? ewaldState_->proposeMove(oldChargedSites, newChargedSites) : 0.0;
    const double deltaU = (uAfter - uBefore) + (electroAfter - electroBefore) + reciprocalDelta;

    const double ratio = gcmc::metropolisRatio(temperatureKelvin_, deltaU);
    if (rng_->nextUniform() < gcmc::clampToProbability(ratio)) {
        if (ewaldState_ != nullptr) {
            ewaldState_->commitMove(oldChargedSites, newChargedSites);
            maybeResyncEwald();
        }
        return;
    }
    for (std::size_t k = 0; k < sites.size(); ++k) {
        setPosition(particles_, sites[k], original[k]);
    }
}

void MonteCarloEngine::attemptInsertion() {
    const auto& matrix = lattice_.matrix();
    const Vec3 a0{matrix[0][0], matrix[0][1], matrix[0][2]};
    const Vec3 a1{matrix[1][0], matrix[1][1], matrix[1][2]};
    const Vec3 a2{matrix[2][0], matrix[2][1], matrix[2][2]};
    const Vec3 referencePoint =
        add(add(scale(a0, rng_->nextUniform()), scale(a1, rng_->nextUniform())),
            scale(a2, rng_->nextUniform()));

    Quaternion orientation{1.0, 0.0, 0.0, 0.0};
    if (species_.sites.size() > 1) {
        orientation =
            randomUnitQuaternion(rng_->nextUniform(), rng_->nextUniform(), rng_->nextUniform());
    }

    std::vector<std::size_t> newSites;
    newSites.reserve(species_.sites.size());
    for (const auto& site : species_.sites) {
        const Vec3 local{site.x, site.y, site.z};
        const Vec3 rotated =
            species_.sites.size() > 1 ? rotateByQuaternion(local, orientation) : local;
        const Vec3 worldPos = add(referencePoint, rotated);
        newSites.push_back(particles_.push_back(worldPos[0], worldPos[1], worldPos[2], 0.0, 0.0,
                                                  0.0, site.mass, site.charge, site.ljSpecies));
    }

    const double dispersionDeltaU =
        moleculeInteractionEnergy(*forceField_, particles_, lattice_, newSites);
    const double electroRealSpace = moleculeElectrostaticRealSpaceEnergy(newSites);
    const auto newChargedSites =
        ewaldState_ != nullptr ? chargedSitesOf(particles_, newSites) : std::vector<forcefield::ChargedSite>{};
    const double reciprocalDelta =
        ewaldState_ != nullptr ? ewaldState_->proposeInsertion(newChargedSites) : 0.0;
    const double deltaU = dispersionDeltaU + electroRealSpace + reciprocalDelta;

    const std::size_t countBefore = molecules_.size();
    const double volume = std::abs(lattice_.volume());
    const double ratio =
        gcmc::insertionRatio(fugacityInternal_, volume, countBefore, temperatureKelvin_, deltaU);

    if (rng_->nextUniform() < gcmc::clampToProbability(ratio)) {
        if (ewaldState_ != nullptr) {
            ewaldState_->commitInsertion(newChargedSites);
        }
        molecules_.push_back(std::move(newSites));
        maybeResyncEwald();
        return;
    }
    // Rejected: the trial sites were just appended at the very end of
    // particles_, so removal is a plain truncation.
    for (std::size_t i = 0; i < species_.sites.size(); ++i) {
        particles_.popBack();
    }
}

void MonteCarloEngine::attemptDeletion() {
    if (molecules_.empty()) {
        return;
    }
    const auto moleculeIndex =
        static_cast<std::size_t>(rng_->nextUniform() * static_cast<double>(molecules_.size()));
    const auto& sites = molecules_[moleculeIndex];

    const double uInteraction = moleculeInteractionEnergy(*forceField_, particles_, lattice_, sites);
    const double electroRealSpace = moleculeElectrostaticRealSpaceEnergy(sites);
    const auto existingChargedSites =
        ewaldState_ != nullptr ? chargedSitesOf(particles_, sites) : std::vector<forcefield::ChargedSite>{};
    const double reciprocalDelta =
        ewaldState_ != nullptr ? ewaldState_->proposeDeletion(existingChargedSites) : 0.0;
    // energy change of REMOVING this molecule: the real-space/dispersion
    // interaction it currently has is lost (negative contribution);
    // reciprocalDelta already carries the correct sign for "this
    // molecule's charges leave" (EwaldIncrementalState::proposeDeletion).
    const double deltaU = -(uInteraction + electroRealSpace) + reciprocalDelta;

    const std::size_t countBefore = molecules_.size();
    const double volume = std::abs(lattice_.volume());
    const double ratio =
        gcmc::deletionRatio(fugacityInternal_, volume, countBefore, temperatureKelvin_, deltaU);

    if (rng_->nextUniform() < gcmc::clampToProbability(ratio)) {
        if (ewaldState_ != nullptr) {
            ewaldState_->commitDeletion(existingChargedSites);
        }
        removeMolecule(moleculeIndex);
        maybeResyncEwald();
    }
}

void MonteCarloEngine::removeMolecule(std::size_t moleculeIndex) {
    std::vector<std::size_t> indices = molecules_[moleculeIndex]; // local copy: fixed up in place
    for (std::size_t k = 0; k < indices.size(); ++k) {
        const std::size_t site = indices[k];
        const std::size_t last = particles_.size() - 1;
        if (site != last) {
            particles_.copyParticle(site, last);
            // Any molecule (including this one's not-yet-removed
            // remaining sites) that referenced index `last` now lives at
            // `site` — fix up every reference before it's read again.
            for (auto& molecule : molecules_) {
                for (auto& s : molecule) {
                    if (s == last) {
                        s = site;
                    }
                }
            }
            for (std::size_t k2 = k + 1; k2 < indices.size(); ++k2) {
                if (indices[k2] == last) {
                    indices[k2] = site;
                }
            }
        }
        particles_.popBack();
    }
    molecules_.erase(molecules_.begin() + static_cast<std::ptrdiff_t>(moleculeIndex));
}

double MonteCarloEngine::widomInsertionHenryCoefficient(std::size_t numTrials) const {
    core::ParticleData framework;
    framework.resize(frameworkCount_);
    for (std::size_t i = 0; i < frameworkCount_; ++i) {
        framework.x[i] = particles_.x[i];
        framework.y[i] = particles_.y[i];
        framework.z[i] = particles_.z[i];
        framework.mass[i] = particles_.mass[i];
        framework.charge[i] = particles_.charge[i];
        framework.species[i] = particles_.species[i];
    }

    const auto& matrix = lattice_.matrix();
    const Vec3 a0{matrix[0][0], matrix[0][1], matrix[0][2]};
    const Vec3 a1{matrix[1][0], matrix[1][1], matrix[1][2]};
    const Vec3 a2{matrix[2][0], matrix[2][1], matrix[2][2]};

    // Widom trials are always hypothetical (never applied to `particles_`),
    // so this uses its own EwaldIncrementalState scoped to the local
    // `framework` copy above -- propose*() only, commit*() never called,
    // exactly matching how every trial's positions are pushed then popped
    // regardless of deltaU.
    std::unique_ptr<forcefield::EwaldIncrementalState> widomEwaldState;
    if (electrostatics_ != nullptr) {
        widomEwaldState = std::make_unique<forcefield::EwaldIncrementalState>(*electrostatics_,
                                                                                framework, lattice_);
    }

    double sumBoltzmann = 0.0;
    for (std::size_t trial = 0; trial < numTrials; ++trial) {
        const Vec3 referencePoint =
            add(add(scale(a0, rng_->nextUniform()), scale(a1, rng_->nextUniform())),
                scale(a2, rng_->nextUniform()));
        Quaternion orientation{1.0, 0.0, 0.0, 0.0};
        if (species_.sites.size() > 1) {
            orientation =
                randomUnitQuaternion(rng_->nextUniform(), rng_->nextUniform(), rng_->nextUniform());
        }

        std::vector<std::size_t> trialSites;
        trialSites.reserve(species_.sites.size());
        for (const auto& site : species_.sites) {
            const Vec3 local{site.x, site.y, site.z};
            const Vec3 rotated =
                species_.sites.size() > 1 ? rotateByQuaternion(local, orientation) : local;
            const Vec3 worldPos = add(referencePoint, rotated);
            trialSites.push_back(framework.push_back(worldPos[0], worldPos[1], worldPos[2], 0.0,
                                                       0.0, 0.0, site.mass, site.charge,
                                                       site.ljSpecies));
        }

        double deltaU = moleculeInteractionEnergy(*forceField_, framework, lattice_, trialSites);
        if (electrostatics_ != nullptr) {
            deltaU +=
                electrostaticRealSpaceInteractionEnergy(*electrostatics_, framework, lattice_,
                                                          trialSites);
            deltaU += widomEwaldState->proposeInsertion(chargedSitesOf(framework, trialSites));
        }
        sumBoltzmann += std::exp(-deltaU / temperatureKelvin_);

        for (std::size_t i = 0; i < species_.sites.size(); ++i) {
            framework.popBack();
        }
    }
    return (sumBoltzmann / static_cast<double>(numTrials)) / temperatureKelvin_;
}

double MonteCarloEngine::electrostaticEnergyDriftForTesting() const {
    if (ewaldState_ == nullptr) {
        return 0.0;
    }
    const double committed = ewaldState_->committedEnergy();
    const double fresh = ewaldState_->recomputeFromScratch(particles_, molecules_);
    return std::abs(fresh - committed);
}

double MonteCarloEngine::electrostaticEnergyForTesting() const {
    return ewaldState_ != nullptr ? ewaldState_->committedEnergy() : 0.0;
}

} // namespace aleator::engines
