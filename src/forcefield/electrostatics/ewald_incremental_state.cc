#include "forcefield/electrostatics/ewald_incremental_state.hpp"

#include <array>
#include <cmath>
#include <numbers>

namespace aleator::forcefield {

namespace {

using Vec3 = std::array<double, 3>;

Vec3 position(const ChargedSite& s) { return {s.x, s.y, s.z}; }

double dot(const Vec3& a, const Vec3& b) { return a[0] * b[0] + a[1] * b[1] + a[2] * b[2]; }
double norm(const Vec3& a) { return std::sqrt(dot(a, a)); }

/// erf(alpha*r)/r with the removable singularity at r -> 0 handled via the
/// analytic limit (same convention as ewald.cc's erfOverR -- duplicated
/// here rather than shared across translation units for a one-line
/// function, matching this codebase's existing precedent of small
/// self-contained per-file math helpers, e.g. monte_carlo_engine.cc's own
/// Vec3/dot/cross rather than reusing ewald.cc's identical-looking ones).
double erfOverR(double r, double alpha) {
    constexpr double kSmallR = 1e-8;
    if (r < kSmallR) {
        return 2.0 * alpha / std::sqrt(std::numbers::pi);
    }
    return std::erf(alpha * r) / r;
}

double pairExclusionSum(const std::vector<ChargedSite>& sites, const core::Lattice& lattice,
                         double alpha) {
    double sum = 0.0;
    for (std::size_t i = 0; i < sites.size(); ++i) {
        for (std::size_t j = i + 1; j < sites.size(); ++j) {
            const double r =
                norm(lattice.minimumImageDisplacement(position(sites[i]), position(sites[j])));
            sum += sites[i].charge * sites[j].charge * erfOverR(r, alpha);
        }
    }
    return sum;
}

double sumChargeSquared(const std::vector<ChargedSite>& sites) {
    double sum = 0.0;
    for (const auto& s : sites) {
        sum += s.charge * s.charge;
    }
    return sum;
}

} // namespace

EwaldIncrementalState::EwaldIncrementalState(const Ewald& ewald,
                                              const core::ParticleData& particles,
                                              const core::Lattice& lattice)
    : ewald_(&ewald),
      lattice_(&lattice),
      terms_(ewald.reciprocalTerms(lattice)),
      volume_(std::abs(lattice.volume())) {
    sRe_.assign(terms_.size(), 0.0);
    sIm_.assign(terms_.size(), 0.0);
    for (std::size_t t = 0; t < terms_.size(); ++t) {
        const auto& k = terms_[t].k;
        double re = 0.0;
        double im = 0.0;
        for (std::size_t i = 0; i < particles.size(); ++i) {
            const double kr =
                k[0] * particles.x[i] + k[1] * particles.y[i] + k[2] * particles.z[i];
            re += particles.charge[i] * std::cos(kr);
            im += particles.charge[i] * std::sin(kr);
        }
        sRe_[t] = re;
        sIm_[t] = im;
    }
    double q2 = 0.0;
    for (double q : particles.charge) {
        q2 += q * q;
    }
    selfEnergy_ = -Ewald::kCoulombConstant * (ewald_->alpha() / std::sqrt(std::numbers::pi)) * q2;
    exclusionEnergy_ = 0.0; // framework carries no intramolecular exclusion, see header doc.
}

double EwaldIncrementalState::committedEnergy() const noexcept {
    double sum = 0.0;
    for (std::size_t t = 0; t < terms_.size(); ++t) {
        sum += terms_[t].weight * (sRe_[t] * sRe_[t] + sIm_[t] * sIm_[t]);
    }
    const double reciprocal =
        sum * Ewald::kCoulombConstant * (2.0 * std::numbers::pi / volume_);
    return reciprocal + selfEnergy_ + exclusionEnergy_;
}

EwaldIncrementalState::Delta EwaldIncrementalState::computeDelta(
    const std::vector<ChargedSite>& removed, const std::vector<ChargedSite>& added) const {
    Delta delta;
    delta.deltaRe.assign(terms_.size(), 0.0);
    delta.deltaIm.assign(terms_.size(), 0.0);

    const double prefactor = Ewald::kCoulombConstant * (2.0 * std::numbers::pi / volume_);
    double weightedDeltaSum = 0.0;
    for (std::size_t t = 0; t < terms_.size(); ++t) {
        const auto& k = terms_[t].k;
        double dRe = 0.0;
        double dIm = 0.0;
        for (const auto& s : added) {
            const double kr = k[0] * s.x + k[1] * s.y + k[2] * s.z;
            dRe += s.charge * std::cos(kr);
            dIm += s.charge * std::sin(kr);
        }
        for (const auto& s : removed) {
            const double kr = k[0] * s.x + k[1] * s.y + k[2] * s.z;
            dRe -= s.charge * std::cos(kr);
            dIm -= s.charge * std::sin(kr);
        }
        const double newRe = sRe_[t] + dRe;
        const double newIm = sIm_[t] + dIm;
        const double oldMagSquared = sRe_[t] * sRe_[t] + sIm_[t] * sIm_[t];
        const double newMagSquared = newRe * newRe + newIm * newIm;
        weightedDeltaSum += terms_[t].weight * (newMagSquared - oldMagSquared);
        delta.deltaRe[t] = dRe;
        delta.deltaIm[t] = dIm;
    }
    delta.reciprocalEnergy = prefactor * weightedDeltaSum;

    const double selfPrefactor =
        -Ewald::kCoulombConstant * (ewald_->alpha() / std::sqrt(std::numbers::pi));
    delta.selfEnergy = selfPrefactor * (sumChargeSquared(added) - sumChargeSquared(removed));

    const double addedExclusion = pairExclusionSum(added, *lattice_, ewald_->alpha());
    const double removedExclusion = pairExclusionSum(removed, *lattice_, ewald_->alpha());
    delta.exclusionEnergy = -Ewald::kCoulombConstant * (addedExclusion - removedExclusion);

    return delta;
}

void EwaldIncrementalState::apply(const Delta& delta) {
    for (std::size_t t = 0; t < terms_.size(); ++t) {
        sRe_[t] += delta.deltaRe[t];
        sIm_[t] += delta.deltaIm[t];
    }
    selfEnergy_ += delta.selfEnergy;
    exclusionEnergy_ += delta.exclusionEnergy;
}

double EwaldIncrementalState::proposeMove(const std::vector<ChargedSite>& oldSites,
                                           const std::vector<ChargedSite>& newSites) const {
    return computeDelta(oldSites, newSites).total();
}

double EwaldIncrementalState::proposeInsertion(const std::vector<ChargedSite>& newSites) const {
    return computeDelta({}, newSites).total();
}

double EwaldIncrementalState::proposeDeletion(
    const std::vector<ChargedSite>& existingSites) const {
    return computeDelta(existingSites, {}).total();
}

void EwaldIncrementalState::commitMove(const std::vector<ChargedSite>& oldSites,
                                        const std::vector<ChargedSite>& newSites) {
    apply(computeDelta(oldSites, newSites));
}

void EwaldIncrementalState::commitInsertion(const std::vector<ChargedSite>& newSites) {
    apply(computeDelta({}, newSites));
}

void EwaldIncrementalState::commitDeletion(const std::vector<ChargedSite>& existingSites) {
    apply(computeDelta(existingSites, {}));
}

double EwaldIncrementalState::recomputeFromScratch(
    const core::ParticleData& particles,
    const std::vector<std::vector<std::size_t>>& molecules) const {
    double reciprocalSum = 0.0;
    for (const auto& term : terms_) {
        double re = 0.0;
        double im = 0.0;
        for (std::size_t i = 0; i < particles.size(); ++i) {
            const double kr = term.k[0] * particles.x[i] + term.k[1] * particles.y[i] +
                               term.k[2] * particles.z[i];
            re += particles.charge[i] * std::cos(kr);
            im += particles.charge[i] * std::sin(kr);
        }
        reciprocalSum += term.weight * (re * re + im * im);
    }
    const double reciprocal =
        reciprocalSum * Ewald::kCoulombConstant * (2.0 * std::numbers::pi / volume_);

    double q2 = 0.0;
    for (double q : particles.charge) {
        q2 += q * q;
    }
    const double self =
        -Ewald::kCoulombConstant * (ewald_->alpha() / std::sqrt(std::numbers::pi)) * q2;

    double exclusion = 0.0;
    for (const auto& molecule : molecules) {
        std::vector<ChargedSite> sites;
        sites.reserve(molecule.size());
        for (std::size_t index : molecule) {
            sites.push_back({particles.x[index], particles.y[index], particles.z[index],
                              particles.charge[index]});
        }
        exclusion +=
            -Ewald::kCoulombConstant * pairExclusionSum(sites, *lattice_, ewald_->alpha());
    }

    return reciprocal + self + exclusion;
}

void EwaldIncrementalState::resync(const core::ParticleData& particles,
                                    const std::vector<std::vector<std::size_t>>& molecules) {
    for (std::size_t t = 0; t < terms_.size(); ++t) {
        const auto& k = terms_[t].k;
        double re = 0.0;
        double im = 0.0;
        for (std::size_t i = 0; i < particles.size(); ++i) {
            const double kr =
                k[0] * particles.x[i] + k[1] * particles.y[i] + k[2] * particles.z[i];
            re += particles.charge[i] * std::cos(kr);
            im += particles.charge[i] * std::sin(kr);
        }
        sRe_[t] = re;
        sIm_[t] = im;
    }

    double q2 = 0.0;
    for (double q : particles.charge) {
        q2 += q * q;
    }
    selfEnergy_ = -Ewald::kCoulombConstant * (ewald_->alpha() / std::sqrt(std::numbers::pi)) * q2;

    double exclusion = 0.0;
    for (const auto& molecule : molecules) {
        std::vector<ChargedSite> sites;
        sites.reserve(molecule.size());
        for (std::size_t index : molecule) {
            sites.push_back({particles.x[index], particles.y[index], particles.z[index],
                              particles.charge[index]});
        }
        exclusion += pairExclusionSum(sites, *lattice_, ewald_->alpha());
    }
    exclusionEnergy_ = -Ewald::kCoulombConstant * exclusion;
}

std::vector<double> EwaldIncrementalState::rawStateForTesting() const {
    std::vector<double> raw;
    raw.reserve(sRe_.size() * 2 + 2);
    raw.insert(raw.end(), sRe_.begin(), sRe_.end());
    raw.insert(raw.end(), sIm_.begin(), sIm_.end());
    raw.push_back(selfEnergy_);
    raw.push_back(exclusionEnergy_);
    return raw;
}

} // namespace aleator::forcefield
