#include "forcefield/pairwise/lennard_jones.hpp"

#include <array>
#include <cmath>
#include <numbers>
#include <stdexcept>
#include <utility>
#include <vector>

namespace aleator::forcefield {

namespace {

/// Raw (untruncated, unshifted) V_LJ(r) = 4*epsilon*[(sigma/r)^12 - (sigma/r)^6].
double rawPotential(double r, const LennardJonesParameters& params) {
    const double sr = params.sigma / r;
    const double sr6 = sr * sr * sr * sr * sr * sr;
    const double sr12 = sr6 * sr6;
    return 4.0 * params.epsilon * (sr12 - sr6);
}

/// dV_LJ/dr = (24*epsilon/r) * [(sigma/r)^6 - 2*(sigma/r)^12].
double rawDerivative(double r, const LennardJonesParameters& params) {
    const double sr = params.sigma / r;
    const double sr6 = sr * sr * sr * sr * sr * sr;
    const double sr12 = sr6 * sr6;
    return (24.0 * params.epsilon / r) * (sr6 - 2.0 * sr12);
}

std::array<double, 3> position(const core::ParticleData& particles, std::size_t i) {
    return {particles.x[i], particles.y[i], particles.z[i]};
}

double norm(const std::array<double, 3>& v) {
    return std::sqrt(v[0] * v[0] + v[1] * v[1] + v[2] * v[2]);
}

} // namespace

LennardJones::LennardJones(std::vector<LennardJonesParameters> perSpecies, double cutoff,
                            LennardJonesTruncation truncation, MixingRule mixingRule)
    : perSpecies_(std::move(perSpecies)), cutoff_(cutoff), truncation_(truncation),
      mixingRule_(mixingRule) {
    if (!(cutoff_ > 0.0)) {
        throw std::invalid_argument("LennardJones: cutoff must be positive");
    }
}

LennardJonesParameters LennardJones::mixedParameters(std::uint32_t speciesI,
                                                       std::uint32_t speciesJ) const {
    const LennardJonesParameters& pi = perSpecies_.at(speciesI);
    if (speciesI == speciesJ) {
        return pi;
    }
    const LennardJonesParameters& pj = perSpecies_.at(speciesJ);

    LennardJonesParameters mixed;
    mixed.epsilon = std::sqrt(pi.epsilon * pj.epsilon); // Berthelot: geometric mean, both rules
    switch (mixingRule_) {
        case MixingRule::LorentzBerthelot:
            mixed.sigma = 0.5 * (pi.sigma + pj.sigma); // Lorentz: arithmetic mean
            break;
        case MixingRule::Geometric:
            mixed.sigma = std::sqrt(pi.sigma * pj.sigma);
            break;
    }
    return mixed;
}

double LennardJones::pairEnergy(double r, const LennardJonesParameters& params) const {
    if (r > cutoff_) {
        return 0.0;
    }
    switch (truncation_) {
        case LennardJonesTruncation::Truncated:
            return rawPotential(r, params);
        case LennardJonesTruncation::Shifted:
            return rawPotential(r, params) - rawPotential(cutoff_, params);
        case LennardJonesTruncation::LinearForceShifted: {
            const double vAtCutoff = rawPotential(cutoff_, params);
            const double dvdrAtCutoff = rawDerivative(cutoff_, params);
            return rawPotential(r, params) - vAtCutoff - (r - cutoff_) * dvdrAtCutoff;
        }
    }
    return 0.0; // unreachable; silences -Wreturn-type on some compilers
}

double LennardJones::pairForceOverR(double r, const LennardJonesParameters& params) const {
    if (r > cutoff_) {
        return 0.0;
    }
    double dvdr = 0.0;
    switch (truncation_) {
        case LennardJonesTruncation::Truncated:
        case LennardJonesTruncation::Shifted:
            // A constant energy shift does not change the derivative.
            dvdr = rawDerivative(r, params);
            break;
        case LennardJonesTruncation::LinearForceShifted:
            dvdr = rawDerivative(r, params) - rawDerivative(cutoff_, params);
            break;
    }
    return -dvdr / r;
}

double LennardJones::computeEnergy(const core::ParticleData& particles, const core::Lattice& lattice,
                                    const core::NeighborList& neighbors) const {
    double energy = 0.0;
    for (const auto& pair : neighbors.pairs()) {
        const auto params =
            mixedParameters(particles.species[pair.i], particles.species[pair.j]);
        const auto d = lattice.minimumImageDisplacement(position(particles, pair.i),
                                                          position(particles, pair.j));
        energy += pairEnergy(norm(d), params);
    }
    return energy;
}

double LennardJones::computeVirial(const core::ParticleData& particles, const core::Lattice& lattice,
                                    const core::NeighborList& neighbors) const {
    double virial = 0.0;
    for (const auto& pair : neighbors.pairs()) {
        const auto params =
            mixedParameters(particles.species[pair.i], particles.species[pair.j]);
        const auto d = lattice.minimumImageDisplacement(position(particles, pair.i),
                                                          position(particles, pair.j));
        const double r = norm(d);
        // W = -sum r * dV/dr; pairForceOverR = (-dV/dr)/r, so
        // -r*dV/dr = pairForceOverR(r) * r^2.
        virial += pairForceOverR(r, params) * r * r;
    }
    return virial;
}

void LennardJones::computeForces(const core::ParticleData& particles, const core::Lattice& lattice,
                                  const core::NeighborList& neighbors, Forces& forcesOut) const {
    for (const auto& pair : neighbors.pairs()) {
        const auto params =
            mixedParameters(particles.species[pair.i], particles.species[pair.j]);
        const auto d = lattice.minimumImageDisplacement(position(particles, pair.i),
                                                          position(particles, pair.j));
        const double r = norm(d);
        const double factor = pairForceOverR(r, params); // (-dV/dr)/r

        // F_i = (dV/dr / r) * d = -factor * d; F_j = -F_i (Newton's third law).
        const double fx = factor * d[0];
        const double fy = factor * d[1];
        const double fz = factor * d[2];
        forcesOut.fx[pair.i] -= fx;
        forcesOut.fy[pair.i] -= fy;
        forcesOut.fz[pair.i] -= fz;
        forcesOut.fx[pair.j] += fx;
        forcesOut.fy[pair.j] += fy;
        forcesOut.fz[pair.j] += fz;
    }
}

double LennardJones::computeParticleEnergy(std::size_t index, const core::ParticleData& particles,
                                            const core::Lattice& lattice,
                                            const std::vector<std::size_t>& excludedIndices) const {
    double energy = 0.0;
    const auto pi = position(particles, index);
    for (std::size_t j = 0; j < particles.size(); ++j) {
        if (j == index) {
            continue;
        }
        bool excluded = false;
        for (std::size_t ex : excludedIndices) {
            if (ex == j) {
                excluded = true;
                break;
            }
        }
        if (excluded) {
            continue;
        }
        const auto params = mixedParameters(particles.species[index], particles.species[j]);
        const auto d = lattice.minimumImageDisplacement(pi, position(particles, j));
        energy += pairEnergy(norm(d), params);
    }
    return energy;
}

namespace {

/// sum_a sum_b N_a * N_b * f(mixedParameters(a, b)), over every ordered
/// pair of species with at least one particle each — the shared structure
/// of tailEnergyCorrection and tailVirialCorrection, which differ only in
/// which per-pair integral f() they use.
template <typename Integral>
double sumOverSpeciesPairs(const LennardJones& lj, const core::ParticleData& particles,
                            Integral&& integral) {
    const auto& species = lj.perSpecies();
    std::vector<std::size_t> counts(species.size(), 0);
    for (auto s : particles.species) {
        counts.at(s) += 1;
    }

    double total = 0.0;
    for (std::size_t a = 0; a < species.size(); ++a) {
        if (counts[a] == 0) {
            continue;
        }
        for (std::size_t b = 0; b < species.size(); ++b) {
            if (counts[b] == 0) {
                continue;
            }
            const auto params =
                lj.mixedParameters(static_cast<std::uint32_t>(a), static_cast<std::uint32_t>(b));
            total += static_cast<double>(counts[a]) * static_cast<double>(counts[b]) *
                     integral(params);
        }
    }
    return total;
}

} // namespace

double LennardJones::tailEnergyCorrection(const core::ParticleData& particles,
                                           const core::Lattice& lattice) const {
    const double volume = std::abs(lattice.volume());
    const double rc = cutoff_;
    const double sum = sumOverSpeciesPairs(*this, particles, [rc](const LennardJonesParameters& p) {
        const double sr3 = p.sigma * p.sigma * p.sigma;
        const double sr9 = sr3 * sr3 * sr3 / (rc * rc * rc * rc * rc * rc * rc * rc * rc);
        const double sr3OverRc3 = sr3 / (rc * rc * rc);
        return p.epsilon * sr3 * ((4.0 / 9.0) * sr9 - (4.0 / 3.0) * sr3OverRc3);
    });
    return (2.0 * std::numbers::pi / volume) * sum;
}

double LennardJones::tailVirialCorrection(const core::ParticleData& particles,
                                           const core::Lattice& lattice) const {
    const double volume = std::abs(lattice.volume());
    const double rc = cutoff_;
    const double sum = sumOverSpeciesPairs(*this, particles, [rc](const LennardJonesParameters& p) {
        const double sr3 = p.sigma * p.sigma * p.sigma;
        const double sr9 = sr3 * sr3 * sr3 / (rc * rc * rc * rc * rc * rc * rc * rc * rc);
        const double sr3OverRc3 = sr3 / (rc * rc * rc);
        return p.epsilon * sr3 * ((16.0 / 3.0) * sr9 - 8.0 * sr3OverRc3);
    });
    return (2.0 * std::numbers::pi / volume) * sum;
}

} // namespace aleator::forcefield
