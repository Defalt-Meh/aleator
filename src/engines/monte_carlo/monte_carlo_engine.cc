#include "engines/monte_carlo/monte_carlo_engine.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <numbers>
#include <numeric>
#include <stdexcept>
#include <utility>

#include "core/exceptions.hpp"
#include "engines/monte_carlo/gcmc_acceptance.hpp"
#include "forcefield/pairwise/lennard_jones.hpp"

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

/// EwaldIncrementalState::resync()/recomputeFromScratch() take plain
/// site-index lists (no species tag -- that class lives in forcefield/,
/// which must not depend on engines/monte_carlo's Molecule type; CLAUDE.md
/// #3's dependency direction is core -> forcefield -> engines, never the
/// reverse). Adapts MonteCarloEngine::Molecule's richer type at the two call
/// sites that need it.
std::vector<std::vector<std::size_t>> siteListsOnly(
    const std::vector<aleator::engines::MonteCarloEngine::Molecule>& molecules) {
    std::vector<std::vector<std::size_t>> siteLists;
    siteLists.reserve(molecules.size());
    for (const auto& molecule : molecules) {
        siteLists.push_back(molecule.sites);
    }
    return siteLists;
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
                                    std::shared_ptr<const forcefield::Ewald> electrostatics,
                                    std::shared_ptr<const FrameworkEnergyGrid> frameworkEnergyGrid)
    : MonteCarloEngine(std::move(frameworkParticles), lattice, std::move(forceField),
                        std::move(rng), temperatureKelvin,
                        std::vector<MoleculeSpecies>{std::move(species)},
                        std::vector<double>{fugacityPascal}, maxTranslationDisplacementAngstrom,
                        maxRotationAngleRadians, std::move(electrostatics),
                        std::move(frameworkEnergyGrid)) {}

MonteCarloEngine::MonteCarloEngine(core::ParticleData frameworkParticles, core::Lattice lattice,
                                    std::shared_ptr<const forcefield::ForceField> forceField,
                                    std::unique_ptr<core::CounterBasedRng> rng,
                                    double temperatureKelvin,
                                    std::vector<MoleculeSpecies> speciesList,
                                    std::vector<double> fugacitiesPascal,
                                    double maxTranslationDisplacementAngstrom,
                                    double maxRotationAngleRadians,
                                    std::shared_ptr<const forcefield::Ewald> electrostatics,
                                    std::shared_ptr<const FrameworkEnergyGrid> frameworkEnergyGrid)
    : ensemble_(Ensemble::Gcmc),
      particles_(std::move(frameworkParticles)),
      lattice_(lattice),
      forceField_(std::move(forceField)),
      rng_(std::move(rng)),
      temperatureKelvin_(temperatureKelvin),
      frameworkCount_(particles_.size()),
      speciesList_(std::move(speciesList)),
      maxTranslation_(maxTranslationDisplacementAngstrom),
      maxRotationAngle_(maxRotationAngleRadians),
      frameworkEnergyGrid_(std::move(frameworkEnergyGrid)),
      electrostatics_(std::move(electrostatics)) {
    if (speciesList_.empty()) {
        throw std::invalid_argument("MonteCarloEngine: GCMC needs at least one species");
    }
    if (fugacitiesPascal.size() != speciesList_.size()) {
        throw std::invalid_argument(
            "MonteCarloEngine: fugacitiesPascal.size() (" + std::to_string(fugacitiesPascal.size()) +
            ") must equal speciesList.size() (" + std::to_string(speciesList_.size()) + ")");
    }
    for (const auto& species : speciesList_) {
        if (species.sites.empty()) {
            throw std::invalid_argument("MonteCarloEngine: every GCMC species must have at least one site");
        }
    }
    if (!(temperatureKelvin_ > 0.0)) {
        throw std::invalid_argument("MonteCarloEngine: temperatureKelvin must be positive");
    }
    fugacitiesInternal_.reserve(fugacitiesPascal.size());
    for (double fugacityPascal : fugacitiesPascal) {
        const double fugacityInternal = fugacityPascal * kPascalToInternal;
        if (!(fugacityInternal > 0.0)) {
            throw std::invalid_argument("MonteCarloEngine: every fugacityPascal must be positive");
        }
        fugacitiesInternal_.push_back(fugacityInternal);
    }

    // CLAUDE.md invariant #10: every trial move below (translation,
    // rotation, insertion, deletion, swap) calls
    // forceField_->computeParticleEnergy() for the moved molecule against
    // the framework. A force field that doesn't support that (e.g. Ewald --
    // see ewald.hpp) must be rejected here, at construction, not discovered
    // as a NotImplemented throw the first time run() is called --
    // "discovering an unsupported combination six hours into a run is a
    // defect."
    if (!forceField_->supportsSingleParticleEnergy()) {
        throw std::invalid_argument(
            "MonteCarloEngine: force field '" + forceField_->name() +
            "' does not support single-particle trial-move energies "
            "(supportsSingleParticleEnergy() == false), which GCMC's "
            "translation/rotation/insertion/deletion/swap moves all require. "
            "This is not yet implemented for this force field -- see "
            "CLAUDE.md section 0.");
    }

    // CLAUDE.md section 5 performance milestone: guest-guest interactions
    // are accelerated via a core::CellList (see moleculeDispersionEnergy()),
    // which needs a fixed cutoff radius to size its grid. That's a narrower
    // requirement than supportsSingleParticleEnergy() alone -- a force field
    // is not required to have a single scalar cutoff() in general (Ewald's
    // real-space term doesn't hard-truncate the way LennardJones does) -- so
    // it's checked here explicitly via the concrete type, rather than added
    // to the ForceField interface for every force field to satisfy.
    // dynamic_cast, not a silent assumption: an unsupported concrete type
    // throws immediately, by name, same as the capability check above
    // (CLAUDE.md invariant #10), instead of being discovered as a crash or
    // wrong answer once run() starts moving molecules.
    const auto* lennardJonesForceField =
        dynamic_cast<const forcefield::LennardJones*>(forceField_.get());
    if (lennardJonesForceField == nullptr) {
        throw std::invalid_argument(
            "MonteCarloEngine: GCMC's guest-guest neighbor acceleration needs a force "
            "field with a single fixed cutoff radius; '" +
            forceField_->name() + "' is not a supported concrete type for this.");
    }
    guestCellListRadius_ = lennardJonesForceField->cutoff();
    frameworkIndices_.resize(frameworkCount_);
    std::iota(frameworkIndices_.begin(), frameworkIndices_.end(), std::size_t{0});
    rebuildGuestCellList();

    if (electrostatics_ != nullptr) {
        // Every adsorbate site's charge must sum to (numerically) zero, PER
        // SPECIES: GCMC inserts/deletes/swaps whole molecules of one
        // species at a time, and if a single molecule carried net charge,
        // every insertion would drift the system's total charge away from
        // neutrality, which Ewald's reciprocal-space k=0 term (silently
        // dropped, see ewald.cc) is only well-defined for. Caught here,
        // once per species, rather than as a slowly-compounding physics
        // error deep into a run.
        constexpr double kNeutralityTolerance = 1e-8;
        for (std::size_t i = 0; i < speciesList_.size(); ++i) {
            double speciesCharge = 0.0;
            for (const auto& site : speciesList_[i].sites) {
                speciesCharge += site.charge;
            }
            if (std::abs(speciesCharge) > kNeutralityTolerance) {
                throw std::invalid_argument(
                    "MonteCarloEngine: GCMC species " + std::to_string(i) +
                    " is not individually charge-neutral (total site charge = " +
                    std::to_string(speciesCharge) +
                    " e); every insertion/deletion/swap would drift the system away from "
                    "neutrality, which electrostatics requires");
            }
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
    const std::vector<std::size_t>& sites, const std::vector<std::size_t>& alsoExclude) const {
    if (electrostatics_ == nullptr) {
        return 0.0;
    }
    if (alsoExclude.empty()) {
        return electrostaticRealSpaceInteractionEnergy(*electrostatics_, particles_, lattice_, sites);
    }
    std::vector<std::size_t> excluded = sites;
    excluded.insert(excluded.end(), alsoExclude.begin(), alsoExclude.end());
    double energy = 0.0;
    for (std::size_t s : sites) {
        energy += electrostatics_->realSpaceParticleEnergy(s, particles_, lattice_, excluded);
    }
    return energy;
}

// CLAUDE.md section 5 performance milestone: replaces the old full-O(N)-scan
// moleculeInteractionEnergy(*forceField_, particles_, lattice_, sites) on
// the GCMC hot path (translation/rotation/insertion/deletion) with:
//  (a) the guest-host (framework) contribution, either an O(1)
//      trilinear-interpolation lookup via frameworkEnergyGrid_ (when
//      given -- see framework_energy_grid.hpp for the accuracy-vs-spacing
//      tradeoff; this part is NOT bit-identical to the direct scan, unlike
//      (b)) or, when frameworkEnergyGrid_ is null, the original direct
//      O(frameworkCount_) scan over frameworkIndices_;
//  (b) the guest-guest contribution, found via guestCellList_ instead of
//      scanning the whole particles_ array.
//
// (b) alone (frameworkEnergyGrid_ == nullptr) is bit-identical to the old
// full scan by construction, not merely close: LennardJones::pairEnergy()
// returns exactly 0.0 for any pair beyond cutoff(), and IEEE754 addition has
// x + 0.0 == x bit-for-bit -- so the old scan's result is unchanged by
// *dropping* any candidate outside cutoff() from the sum, as long as every
// candidate that DOES contribute a nonzero term is visited in the same
// relative order as the old ascending (j = 0, 1, ..., N-1) scan visited it
// (floating-point addition is not associative, so reordering nonzero terms,
// not just omitting zero ones, can change the last bit). frameworkIndices_
// is already in ascending order; guestCandidates is explicitly sorted back
// into ascending order below (CellList's bin traversal order is not index
// order) precisely to preserve that. See
// tests/validation/test_gcmc_cell_list_bit_identical.cc.
double MonteCarloEngine::moleculeDispersionEnergy(const std::vector<std::size_t>& sites,
                                                   const std::vector<std::size_t>& alsoExclude) const {
    std::vector<std::size_t> excluded;
    const std::vector<std::size_t>* excludedPtr = &sites;
    if (!alsoExclude.empty()) {
        excluded = sites;
        excluded.insert(excluded.end(), alsoExclude.begin(), alsoExclude.end());
        excludedPtr = &excluded;
    }

    double energy = 0.0;
    for (std::size_t s : sites) {
        if (frameworkEnergyGrid_ != nullptr) {
            energy += frameworkEnergyGrid_->interpolate(
                particles_.species[s], {particles_.x[s], particles_.y[s], particles_.z[s]});
        }

        std::vector<std::size_t> guestCandidates;
        guestCellList_.forEachIndexNear(
            lattice_, {particles_.x[s], particles_.y[s], particles_.z[s]},
            [&](std::size_t j) { guestCandidates.push_back(j); });
        std::sort(guestCandidates.begin(), guestCandidates.end());

        std::vector<std::size_t> candidates;
        if (frameworkEnergyGrid_ == nullptr) {
            candidates = frameworkIndices_;
        }
        candidates.insert(candidates.end(), guestCandidates.begin(), guestCandidates.end());

        energy += forceField_->computeParticleEnergyOverCandidates(s, particles_, lattice_,
                                                                     candidates, *excludedPtr);
    }
    return energy;
}

void MonteCarloEngine::rebuildGuestCellList() {
    guestCellList_.build(particles_, lattice_, guestCellListRadius_, frameworkCount_);
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
        ewaldState_->resync(particles_, siteListsOnly(molecules_));
        movesSinceEwaldResync_ = 0;
    }
}

void MonteCarloEngine::run(std::size_t numSteps) {
    if (ensemble_ != Ensemble::Gcmc) {
        throw aleator::NotImplemented("MonteCarloEngine::run (NVT/NPT)");
    }
    // Single-species: EXACT pre-mixture 4-way split (no swap branch even
    // present in the code path) -- see class doc comment on bit-identical
    // backward compatibility. Multi-species: 5-way split including swap.
    const bool mixture = speciesList_.size() >= 2;
    for (std::size_t step = 0; step < numSteps; ++step) {
        const double pick = rng_->nextUniform();
        if (!mixture) {
            if (pick < 0.25) {
                attemptTranslation();
            } else if (pick < 0.5) {
                attemptRotation();
            } else if (pick < 0.75) {
                attemptInsertion();
            } else {
                attemptDeletion();
            }
        } else {
            if (pick < 0.2) {
                attemptTranslation();
            } else if (pick < 0.4) {
                attemptRotation();
            } else if (pick < 0.6) {
                attemptInsertion();
            } else if (pick < 0.8) {
                attemptDeletion();
            } else {
                attemptSwap();
            }
        }
    }
}

void MonteCarloEngine::attemptTranslation() {
    if (molecules_.empty()) {
        return;
    }
    const auto moleculeIndex =
        static_cast<std::size_t>(rng_->nextUniform() * static_cast<double>(molecules_.size()));
    const auto& sites = molecules_[moleculeIndex].sites;

    const double uBefore = moleculeDispersionEnergy(sites);
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

    const double uAfter = moleculeDispersionEnergy(sites);
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
        rebuildGuestCellList();
        return; // accepted, keep the new positions
    }
    for (std::size_t k = 0; k < sites.size(); ++k) {
        setPosition(particles_, sites[k], original[k]);
    }
}

void MonteCarloEngine::attemptRotation() {
    if (molecules_.empty()) {
        return;
    }
    // Single-species, single-site: EXACT pre-mixture short-circuit (no
    // molecule-index draw at all) -- see class doc comment. A mixture can't
    // know whether rotation is a no-op without first knowing which
    // molecule got picked, since different species can have different site
    // counts; that's the branch below.
    if (speciesList_.size() == 1 && speciesList_[0].sites.size() < 2) {
        return;
    }
    const auto moleculeIndex =
        static_cast<std::size_t>(rng_->nextUniform() * static_cast<double>(molecules_.size()));
    const Molecule& molecule = molecules_[moleculeIndex];
    const MoleculeSpecies& thisSpecies = speciesList_[molecule.speciesIndex];
    if (thisSpecies.sites.size() < 2) {
        return; // this particular molecule's species is single-site; nothing to rotate
    }
    const auto& sites = molecule.sites;

    const double uBefore = moleculeDispersionEnergy(sites);
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

    const double uAfter = moleculeDispersionEnergy(sites);
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
        rebuildGuestCellList();
        return;
    }
    for (std::size_t k = 0; k < sites.size(); ++k) {
        setPosition(particles_, sites[k], original[k]);
    }
}

void MonteCarloEngine::attemptInsertion() {
    // Single-species: EXACT pre-mixture behavior -- speciesIndex is always
    // 0, and (critically) NO RNG draw is consumed choosing it, so every
    // subsequent draw (referencePoint, orientation) lands in the same
    // position in the stream as before this milestone. See class doc
    // comment.
    const std::size_t speciesIndex =
        speciesList_.size() == 1
            ? std::size_t{0}
            : static_cast<std::size_t>(rng_->nextUniform() * static_cast<double>(speciesList_.size()));
    const MoleculeSpecies& species = speciesList_[speciesIndex];

    const auto& matrix = lattice_.matrix();
    const Vec3 a0{matrix[0][0], matrix[0][1], matrix[0][2]};
    const Vec3 a1{matrix[1][0], matrix[1][1], matrix[1][2]};
    const Vec3 a2{matrix[2][0], matrix[2][1], matrix[2][2]};
    const Vec3 referencePoint =
        add(add(scale(a0, rng_->nextUniform()), scale(a1, rng_->nextUniform())),
            scale(a2, rng_->nextUniform()));

    Quaternion orientation{1.0, 0.0, 0.0, 0.0};
    if (species.sites.size() > 1) {
        orientation =
            randomUnitQuaternion(rng_->nextUniform(), rng_->nextUniform(), rng_->nextUniform());
    }

    std::vector<std::size_t> newSites;
    newSites.reserve(species.sites.size());
    for (const auto& site : species.sites) {
        const Vec3 local{site.x, site.y, site.z};
        const Vec3 rotated = species.sites.size() > 1 ? rotateByQuaternion(local, orientation) : local;
        const Vec3 worldPos = add(referencePoint, rotated);
        newSites.push_back(particles_.push_back(worldPos[0], worldPos[1], worldPos[2], 0.0, 0.0,
                                                  0.0, site.mass, site.charge, site.ljSpecies));
    }

    const double dispersionDeltaU = moleculeDispersionEnergy(newSites);
    const double electroRealSpace = moleculeElectrostaticRealSpaceEnergy(newSites);
    const auto newChargedSites =
        ewaldState_ != nullptr ? chargedSitesOf(particles_, newSites) : std::vector<forcefield::ChargedSite>{};
    const double reciprocalDelta =
        ewaldState_ != nullptr ? ewaldState_->proposeInsertion(newChargedSites) : 0.0;
    const double deltaU = dispersionDeltaU + electroRealSpace + reciprocalDelta;

    const std::size_t countBefore = moleculeCountOfSpecies(speciesIndex);
    const double volume = std::abs(lattice_.volume());
    const double ratio = gcmc::insertionRatio(fugacitiesInternal_[speciesIndex], volume, countBefore,
                                               temperatureKelvin_, deltaU);

    if (rng_->nextUniform() < gcmc::clampToProbability(ratio)) {
        if (ewaldState_ != nullptr) {
            ewaldState_->commitInsertion(newChargedSites);
        }
        molecules_.push_back(Molecule{static_cast<std::uint32_t>(speciesIndex), std::move(newSites)});
        maybeResyncEwald();
        rebuildGuestCellList();
        return;
    }
    // Rejected: the trial sites were just appended at the very end of
    // particles_, so removal is a plain truncation.
    for (std::size_t i = 0; i < species.sites.size(); ++i) {
        particles_.popBack();
    }
}

void MonteCarloEngine::attemptDeletion() {
    if (molecules_.empty()) {
        return;
    }
    std::size_t moleculeIndex;
    std::size_t speciesIndex;
    std::size_t countBefore;
    if (speciesList_.size() == 1) {
        // EXACT pre-mixture behavior/draw sequence -- see class doc comment.
        moleculeIndex =
            static_cast<std::size_t>(rng_->nextUniform() * static_cast<double>(molecules_.size()));
        speciesIndex = 0;
        countBefore = molecules_.size();
    } else {
        speciesIndex = static_cast<std::size_t>(rng_->nextUniform() *
                                                  static_cast<double>(speciesList_.size()));
        const auto indices = moleculeIndicesOfSpecies(speciesIndex);
        if (indices.empty()) {
            return; // nothing of this species to delete; a real, deterministic no-op
        }
        moleculeIndex =
            indices[static_cast<std::size_t>(rng_->nextUniform() * static_cast<double>(indices.size()))];
        countBefore = indices.size();
    }
    const auto& sites = molecules_[moleculeIndex].sites;

    const double uInteraction = moleculeDispersionEnergy(sites);
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

    const double volume = std::abs(lattice_.volume());
    const double ratio = gcmc::deletionRatio(fugacitiesInternal_[speciesIndex], volume, countBefore,
                                              temperatureKelvin_, deltaU);

    if (rng_->nextUniform() < gcmc::clampToProbability(ratio)) {
        if (ewaldState_ != nullptr) {
            ewaldState_->commitDeletion(existingChargedSites);
        }
        removeMolecule(moleculeIndex);
        maybeResyncEwald();
        rebuildGuestCellList();
    }
}

void MonteCarloEngine::attemptSwap() {
    if (molecules_.empty() || speciesList_.size() < 2) {
        return;
    }
    const auto fromSpeciesIndex =
        static_cast<std::size_t>(rng_->nextUniform() * static_cast<double>(speciesList_.size()));
    const auto fromIndices = moleculeIndicesOfSpecies(fromSpeciesIndex);
    if (fromIndices.empty()) {
        return; // nothing of this species present to convert; a real, deterministic no-op
    }
    const auto moleculeIndex =
        fromIndices[static_cast<std::size_t>(rng_->nextUniform() *
                                              static_cast<double>(fromIndices.size()))];
    // Target species uniformly among the OTHER (speciesList_.size()-1)
    // species: draw uniformly over [0, S-1) then skip fromSpeciesIndex --
    // maps onto {0,...,S-1} \ {fromSpeciesIndex} uniformly. Matches
    // swapRatio's derivation, which assumes species selection is symmetric
    // between the forward (A->B) and reverse (B->A) moves.
    auto toSpeciesIndex = static_cast<std::size_t>(
        rng_->nextUniform() * static_cast<double>(speciesList_.size() - 1));
    if (toSpeciesIndex >= fromSpeciesIndex) {
        ++toSpeciesIndex;
    }

    // Captured BEFORE any structural change: molecules_[moleculeIndex] gets
    // reassigned to the new species below, so `oldSites` (the vanishing
    // molecule's current site indices) must be a value copy, not a
    // reference into molecules_ that reassignment would invalidate.
    const std::vector<std::size_t> oldSites = molecules_[moleculeIndex].sites;

    const double uBefore = moleculeDispersionEnergy(oldSites);
    const double electroBefore = moleculeElectrostaticRealSpaceEnergy(oldSites);
    const auto oldChargedSites =
        ewaldState_ != nullptr ? chargedSitesOf(particles_, oldSites) : std::vector<forcefield::ChargedSite>{};

    // New species placed at the OLD molecule's centroid, fresh random
    // orientation -- physically reasonable (higher acceptance than a fresh
    // uniform-random position) and, per swapRatio's derivation, the
    // acceptance-ratio PREFACTOR does not depend on where the new molecule
    // is placed (only deltaU does) -- see gcmc_acceptance.hpp.
    const Vec3 referencePoint = moleculeCentroid(particles_, oldSites);
    const MoleculeSpecies& toSpecies = speciesList_[toSpeciesIndex];
    Quaternion orientation{1.0, 0.0, 0.0, 0.0};
    if (toSpecies.sites.size() > 1) {
        orientation =
            randomUnitQuaternion(rng_->nextUniform(), rng_->nextUniform(), rng_->nextUniform());
    }
    std::vector<std::size_t> newSites;
    newSites.reserve(toSpecies.sites.size());
    for (const auto& site : toSpecies.sites) {
        const Vec3 local{site.x, site.y, site.z};
        const Vec3 rotated =
            toSpecies.sites.size() > 1 ? rotateByQuaternion(local, orientation) : local;
        const Vec3 worldPos = add(referencePoint, rotated);
        newSites.push_back(particles_.push_back(worldPos[0], worldPos[1], worldPos[2], 0.0, 0.0,
                                                  0.0, site.mass, site.charge, site.ljSpecies));
    }

    // The old molecule's sites are still physically present in particles_
    // at this point (not yet removed) -- moleculeDispersionEnergy/
    // moleculeElectrostaticRealSpaceEnergy's `alsoExclude` parameter keeps
    // them from spuriously contributing to the NEW species' trial energy
    // (they occupy, or are very close to, the same position, so without
    // this exclusion the trial would see a large, physically meaningless
    // self-interaction that will not exist in the final, post-swap state).
    const double uAfter = moleculeDispersionEnergy(newSites, oldSites);
    const double electroAfter = moleculeElectrostaticRealSpaceEnergy(newSites, oldSites);
    const auto newChargedSites =
        ewaldState_ != nullptr ? chargedSitesOf(particles_, newSites) : std::vector<forcefield::ChargedSite>{};
    // A swap is exactly EwaldIncrementalState's general "remove these
    // charges, add these charges" delta (computeDelta's doc comment: "an
    // insertion is (removed={}, added=newSites); a deletion is
    // (removed=existingSites, added={})" -- a swap is simply both at once,
    // proposeMove/commitMove's own general form, not a special case of it).
    const double reciprocalDelta =
        ewaldState_ != nullptr ? ewaldState_->proposeMove(oldChargedSites, newChargedSites) : 0.0;
    const double deltaU = (uAfter - uBefore) + (electroAfter - electroBefore) + reciprocalDelta;

    const std::size_t countFromBefore = fromIndices.size();
    const std::size_t countToBefore = moleculeCountOfSpecies(toSpeciesIndex);
    const double ratio =
        gcmc::swapRatio(fugacitiesInternal_[fromSpeciesIndex], fugacitiesInternal_[toSpeciesIndex],
                         countFromBefore, countToBefore, temperatureKelvin_, deltaU);

    if (rng_->nextUniform() < gcmc::clampToProbability(ratio)) {
        if (ewaldState_ != nullptr) {
            ewaldState_->commitMove(oldChargedSites, newChargedSites);
        }
        // Reassign to the new species FIRST, so removeParticleIndices'
        // generic fix-up loop (which scans molecules_) protects the NEW
        // sites' indices if removing the OLD sites happens to relocate any
        // of them via swap-with-last -- the same already-validated
        // mechanism removeMolecule() relies on, not a special case.
        //
        // CRITICAL ORDERING BUG this fixes (caught by
        // tests/validation/test_gcmc_charged_mixture_drift.cc, a real
        // ~5% relative drift explosion, not gradual accumulation):
        // maybeResyncEwald() reads particles_/molecules_ DIRECTLY
        // (EwaldIncrementalState::resync()), so it must never run while
        // particles_ transiently holds BOTH the old (about-to-be-removed)
        // AND new (just-appended) molecule's sites at once, or while
        // molecules_[moleculeIndex] still points at the old sites -- a
        // resync landing in that window would recompute the cache against
        // a physically inconsistent extra-molecule state. Reassigning
        // molecules_ and physically removing the old sites BEFORE calling
        // maybeResyncEwald() (moved below, past both) closes that window.
        molecules_[moleculeIndex] =
            Molecule{static_cast<std::uint32_t>(toSpeciesIndex), std::move(newSites)};
        removeParticleIndices(oldSites);
        maybeResyncEwald();
        rebuildGuestCellList();
        return;
    }
    // Rejected: the trial (new-species) sites were just appended at the
    // very end of particles_; the old molecule is untouched. Plain
    // truncation, same as attemptInsertion()'s rejection path.
    for (std::size_t i = 0; i < toSpecies.sites.size(); ++i) {
        particles_.popBack();
    }
}

std::size_t MonteCarloEngine::moleculeCountOfSpecies(std::size_t speciesIndex) const {
    std::size_t count = 0;
    for (const auto& molecule : molecules_) {
        if (molecule.speciesIndex == speciesIndex) {
            ++count;
        }
    }
    return count;
}

std::vector<std::size_t> MonteCarloEngine::moleculeIndicesOfSpecies(std::size_t speciesIndex) const {
    std::vector<std::size_t> indices;
    for (std::size_t j = 0; j < molecules_.size(); ++j) {
        if (molecules_[j].speciesIndex == speciesIndex) {
            indices.push_back(j);
        }
    }
    return indices;
}

void MonteCarloEngine::removeMolecule(std::size_t moleculeIndex) {
    // Local copy of the sites to remove, taken before erasing the
    // molecules_ entry (removeParticleIndices' own fix-up loop scans
    // molecules_, so erasing first means it never wastes work fixing up an
    // entry that's about to disappear anyway -- harmless either order, this
    // one's just slightly less wasteful).
    std::vector<std::size_t> indices = molecules_[moleculeIndex].sites;
    molecules_.erase(molecules_.begin() + static_cast<std::ptrdiff_t>(moleculeIndex));
    removeParticleIndices(std::move(indices));
}

void MonteCarloEngine::removeParticleIndices(std::vector<std::size_t> indices) {
    for (std::size_t k = 0; k < indices.size(); ++k) {
        const std::size_t site = indices[k];
        const std::size_t last = particles_.size() - 1;
        if (site != last) {
            particles_.copyParticle(site, last);
            // Any molecule that referenced index `last` now lives at
            // `site` — fix up every reference before it's read again.
            for (auto& molecule : molecules_) {
                for (auto& s : molecule.sites) {
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
}

double MonteCarloEngine::widomInsertionHenryCoefficient(std::size_t numTrials,
                                                         std::size_t speciesIndex) const {
    const MoleculeSpecies& species = speciesList_[speciesIndex];
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
        if (species.sites.size() > 1) {
            orientation =
                randomUnitQuaternion(rng_->nextUniform(), rng_->nextUniform(), rng_->nextUniform());
        }

        std::vector<std::size_t> trialSites;
        trialSites.reserve(species.sites.size());
        for (const auto& site : species.sites) {
            const Vec3 local{site.x, site.y, site.z};
            const Vec3 rotated =
                species.sites.size() > 1 ? rotateByQuaternion(local, orientation) : local;
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

        for (std::size_t i = 0; i < species.sites.size(); ++i) {
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
    const double fresh = ewaldState_->recomputeFromScratch(particles_, siteListsOnly(molecules_));
    return std::abs(fresh - committed);
}

double MonteCarloEngine::electrostaticEnergyForTesting() const {
    return ewaldState_ != nullptr ? ewaldState_->committedEnergy() : 0.0;
}

} // namespace aleator::engines
