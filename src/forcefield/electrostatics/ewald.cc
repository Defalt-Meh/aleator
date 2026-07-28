#include "forcefield/electrostatics/ewald.hpp"

#include <array>
#include <cmath>
#include <numbers>
#include <stdexcept>
#include <string>
#include <unordered_set>
#include <utility>

namespace aleator::forcefield {

namespace {

using Vec3 = std::array<double, 3>;

Vec3 add(const Vec3& a, const Vec3& b) { return {a[0] + b[0], a[1] + b[1], a[2] + b[2]}; }
Vec3 scale(const Vec3& a, double s) { return {a[0] * s, a[1] * s, a[2] * s}; }
double dot(const Vec3& a, const Vec3& b) { return a[0] * b[0] + a[1] * b[1] + a[2] * b[2]; }
Vec3 cross(const Vec3& a, const Vec3& b) {
    return {a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0]};
}
double norm(const Vec3& a) { return std::sqrt(dot(a, a)); }

Vec3 position(const core::ParticleData& particles, std::size_t i) {
    return {particles.x[i], particles.y[i], particles.z[i]};
}

/// Reciprocal lattice vectors b_i = 2*pi*(a_j x a_k)/V, satisfying
/// a_i . b_j = 2*pi*delta_ij. Not exposed by core::Lattice (it only needs
/// real-space geometry for everything else in the codebase so far), so
/// Ewald computes them itself from the real-space matrix.
std::array<Vec3, 3> reciprocalVectors(const core::Lattice& lattice) {
    const auto& m = lattice.matrix();
    const Vec3 a0{m[0][0], m[0][1], m[0][2]};
    const Vec3 a1{m[1][0], m[1][1], m[1][2]};
    const Vec3 a2{m[2][0], m[2][1], m[2][2]};
    const double volume = lattice.volume(); // signed; consistent with the cross products below
    const double twoPiOverV = 2.0 * std::numbers::pi / volume;
    return {
        scale(cross(a1, a2), twoPiOverV),
        scale(cross(a2, a0), twoPiOverV),
        scale(cross(a0, a1), twoPiOverV),
    };
}

/// Real-space part of a unit-charge-pair interaction, radial factor such
/// that the pair energy is q_i*q_j*realSpaceKernel(r, alpha):
///     erfc(alpha*r) / r
/// This is the short-ranged piece of 1/r = erf(alpha*r)/r + erfc(alpha*r)/r
/// — summed only over neighbors within realSpaceCutoff_, using the
/// minimum-image convention (valid because cutoff <= L_perp/2 is enforced
/// wherever the neighbor list was built).
double realSpaceKernel(double r, double alpha) { return std::erfc(alpha * r) / r; }

/// (-d/dr)[erfc(alpha*r)/r] / r — the "force over r" factor for the
/// real-space kernel, i.e. multiplying by the (i -> j) minimum-image
/// displacement vector gives the force on i from j directly (same
/// convention as LennardJones::pairForceOverR).
///     d/dr[erfc(a*r)/r] = -erfc(a*r)/r^2 - (2*a/sqrt(pi)) * exp(-a^2*r^2) / r
/// so (-d/dr)/r = erfc(a*r)/r^3 + (2*a/sqrt(pi)) * exp(-a^2*r^2) / r^2.
double realSpaceForceOverR(double r, double alpha) {
    const double r2 = r * r;
    return std::erfc(alpha * r) / (r2 * r) +
           (2.0 * alpha / std::sqrt(std::numbers::pi)) * std::exp(-alpha * alpha * r2) / r2;
}

/// erf(alpha*r)/r, with the removable singularity at r -> 0 handled via
/// the analytic limit erf(x)/x -> 2/sqrt(pi) as x -> 0, i.e.
/// erf(alpha*r)/r -> 2*alpha/sqrt(pi). Used by the exclusion correction;
/// bonded pairs are never literally coincident in practice, but the limit
/// costs nothing to handle correctly.
double erfOverR(double r, double alpha) {
    constexpr double kSmallR = 1e-8;
    if (r < kSmallR) {
        return 2.0 * alpha / std::sqrt(std::numbers::pi);
    }
    return std::erf(alpha * r) / r;
}

/// (-d/dr)[erf(alpha*r)/r] / r, the exclusion-correction "force over r"
/// factor (see ewald.hpp: the exclusion energy term subtracts
/// q_i*q_j*erf(alpha*r)/r for excluded pairs, so its force contribution is
/// the negative gradient of that subtracted term).
///     d/dr[erf(a*r)/r] = (2*a/sqrt(pi)) * exp(-a^2*r^2) / r - erf(a*r)/r^2
/// so (-d/dr)/r = erf(a*r)/r^3 - (2*a/sqrt(pi)) * exp(-a^2*r^2) / r^2.
/// At r -> 0: erf(a*r)/r^3 -> (2*a/sqrt(pi))*(1/r^2) - (2*a^3)/(3*sqrt(pi))
/// (from the Taylor series erf(x) = 2/sqrt(pi)*(x - x^3/3 + ...)), which
/// exactly cancels the -(2*a/sqrt(pi))*exp(...)/r^2 term's leading 1/r^2
/// piece, leaving a finite limit; excluded pairs are bonded atoms and
/// never coincident in any real input, so this limit is not handled
/// separately here (unlike erfOverR above, which the exclusion *energy*
/// does need).
double exclusionForceOverR(double r, double alpha) {
    const double r2 = r * r;
    return std::erf(alpha * r) / (r2 * r) -
           (2.0 * alpha / std::sqrt(std::numbers::pi)) * std::exp(-alpha * alpha * r2) / r2;
}

struct ExclusionSet {
    std::unordered_set<std::size_t> encoded;

    static std::size_t encode(std::size_t i, std::size_t j) {
        if (i > j) {
            std::swap(i, j);
        }
        // Fine as long as no index exceeds 2^32 particles.
        return (i << 32) | j;
    }

    explicit ExclusionSet(const std::vector<std::pair<std::size_t, std::size_t>>& pairs) {
        encoded.reserve(pairs.size());
        for (const auto& [i, j] : pairs) {
            encoded.insert(encode(i, j));
        }
    }

    [[nodiscard]] bool contains(std::size_t i, std::size_t j) const {
        return encoded.contains(encode(i, j));
    }
};

/// One reciprocal-lattice vector's contribution, precomputed once and
/// shared between computeEnergy() and computeForces() so both read from
/// the exact same structure-factor evaluation.
struct ReciprocalTerm {
    Vec3 k;
    double weight; // exp(-k^2 / (4*alpha^2)) / k^2
};

/// Reciprocal-lattice index vectors are truncated to a SPHERE in index
/// space (n0^2 + n1^2 + n2^2 < kMax^2 + 2), not a cube (|n_i| <= kMax
/// independently) — a cubic truncation is anisotropic, pulling in "corner"
/// vectors (e.g. (kMax,kMax,kMax), with |n|^2 = 3*kMax^2) that are much
/// further from the origin than vectors near the axes at the same kMax,
/// while omitting comparably-distant vectors that a spherical cutoff would
/// include. The "+2" and the STRICT inequality match the NIST SRSW SPC/E
/// reference calculations' convention exactly (see
/// tests/validation/test_ewald_spce_water.cc); confirmed empirically —
/// using <= here (including the n^2 == kMax^2+2 shell) reproduced NIST's
/// E_real/E_self/E_intra to ~1e-6 but E_fourier consistently ~4-6e-4 too
/// large across all four reference configurations, because for kMax=5 the
/// n^2 == 27 shell alone is 32 extra vectors: the 24 permutations of
/// (+-5,+-1,+-1) and 8 of (+-3,+-3,+-3). Switching to strict < resolved it.
std::vector<ReciprocalTerm> buildReciprocalTerms(const core::Lattice& lattice, double alpha,
                                                  int kMax) {
    const auto [b0, b1, b2] = reciprocalVectors(lattice);
    const double fourAlphaSquared = 4.0 * alpha * alpha;
    const int kMaxSquaredBound = kMax * kMax + 2;

    std::vector<ReciprocalTerm> terms;
    terms.reserve(static_cast<std::size_t>((2 * kMax + 1) * (2 * kMax + 1) * (2 * kMax + 1)));
    for (int n0 = -kMax; n0 <= kMax; ++n0) {
        for (int n1 = -kMax; n1 <= kMax; ++n1) {
            for (int n2 = -kMax; n2 <= kMax; ++n2) {
                if (n0 == 0 && n1 == 0 && n2 == 0) {
                    continue; // k=0 term is the (dropped) neutralizing background
                }
                if (n0 * n0 + n1 * n1 + n2 * n2 >= kMaxSquaredBound) {
                    continue;
                }
                const Vec3 k = add(add(scale(b0, n0), scale(b1, n1)), scale(b2, n2));
                const double k2 = dot(k, k);
                terms.push_back({k, std::exp(-k2 / fourAlphaSquared) / k2});
            }
        }
    }
    return terms;
}

/// Structure factor S(k) = sum_i q_i * exp(i k.r_i), split into real and
/// imaginary parts (charges are real, so this is just Re/Im accumulation,
/// no complex<double> needed).
struct StructureFactor {
    double re = 0.0;
    double im = 0.0;
};

StructureFactor structureFactor(const core::ParticleData& particles, const Vec3& k) {
    StructureFactor s;
    for (std::size_t i = 0; i < particles.size(); ++i) {
        const double kr = dot(k, position(particles, i));
        s.re += particles.charge[i] * std::cos(kr);
        s.im += particles.charge[i] * std::sin(kr);
    }
    return s;
}

void validateChargeNeutrality(const core::ParticleData& particles) {
    double totalCharge = 0.0;
    for (double q : particles.charge) {
        totalCharge += q;
    }
    // The k=0 reciprocal term (dropped above) is only finite for a neutral
    // system; a nonzero net charge needs an explicit neutralizing
    // background to even be well-defined, which is not implemented.
    constexpr double kNeutralityTolerance = 1e-8;
    if (std::abs(totalCharge) > kNeutralityTolerance) {
        throw std::invalid_argument(
            "Ewald: system is not charge-neutral (total charge = " + std::to_string(totalCharge) +
            " e); the reciprocal-space k=0 term is only well-defined for a neutral system");
    }
}

} // namespace

Ewald::Ewald(double alpha, double realSpaceCutoff, int kMax,
             std::vector<std::pair<std::size_t, std::size_t>> exclusions)
    : alpha_(alpha), realSpaceCutoff_(realSpaceCutoff), kMax_(kMax),
      exclusions_(std::move(exclusions)) {
    if (!(alpha_ > 0.0)) {
        throw std::invalid_argument("Ewald: alpha must be positive");
    }
    if (!(realSpaceCutoff_ > 0.0)) {
        throw std::invalid_argument("Ewald: realSpaceCutoff must be positive");
    }
    if (kMax_ < 1) {
        throw std::invalid_argument("Ewald: kMax must be at least 1");
    }
}

double Ewald::computeEnergy(const core::ParticleData& particles, const core::Lattice& lattice,
                             const core::NeighborList& neighbors) const {
    return computeEnergyBreakdown(particles, lattice, neighbors).total();
}

EwaldEnergyBreakdown Ewald::computeEnergyBreakdown(const core::ParticleData& particles,
                                                    const core::Lattice& lattice,
                                                    const core::NeighborList& neighbors) const {
    validateChargeNeutrality(particles);
    const ExclusionSet excluded(exclusions_);
    EwaldEnergyBreakdown breakdown;

    // --- U_real: short-ranged erfc(alpha*r)/r sum over neighbor pairs,
    // minimum image, skipping excluded (e.g. intramolecular) pairs. ---
    double uReal = 0.0;
    for (const auto& pair : neighbors.pairs()) {
        if (excluded.contains(pair.i, pair.j)) {
            continue;
        }
        const double r = norm(lattice.minimumImageDisplacement(position(particles, pair.i),
                                                                 position(particles, pair.j)));
        if (r > realSpaceCutoff_) {
            continue; // neighbor list may have been built at a larger cutoff
        }
        uReal += particles.charge[pair.i] * particles.charge[pair.j] * realSpaceKernel(r, alpha_);
    }
    breakdown.real = uReal * kCoulombConstant;

    // --- U_reciprocal: (2*pi/V) * sum_{k!=0} exp(-k^2/4a^2)/k^2 * |S(k)|^2.
    // ---
    const double volume = std::abs(lattice.volume());
    const auto reciprocalTerms = buildReciprocalTerms(lattice, alpha_, kMax_);
    double uReciprocal = 0.0;
    for (const auto& term : reciprocalTerms) {
        const auto s = structureFactor(particles, term.k);
        uReciprocal += term.weight * (s.re * s.re + s.im * s.im);
    }
    breakdown.reciprocal = uReciprocal * kCoulombConstant * (2.0 * std::numbers::pi / volume);

    // --- U_self: removes each charge's spurious interaction with its own
    // neutralizing Gaussian screening cloud, implicitly present in
    // U_reciprocal. ---
    double sumChargeSquared = 0.0;
    for (double q : particles.charge) {
        sumChargeSquared += q * q;
    }
    breakdown.self =
        -kCoulombConstant * (alpha_ / std::sqrt(std::numbers::pi)) * sumChargeSquared;

    // --- U_exclusion: removes the erf(alpha*r)/r contribution that
    // U_reciprocal implicitly includes for excluded (e.g. intramolecular)
    // pairs, which U_real never added in the first place (they were
    // skipped above). ---
    double uExclusion = 0.0;
    for (const auto& [i, j] : exclusions_) {
        const double r =
            norm(lattice.minimumImageDisplacement(position(particles, i), position(particles, j)));
        uExclusion += particles.charge[i] * particles.charge[j] * erfOverR(r, alpha_);
    }
    breakdown.exclusion = -uExclusion * kCoulombConstant;

    // --- U_surface: exactly 0 under tinfoil (conducting) boundary
    // conditions — see class doc comment in ewald.hpp. Written out
    // explicitly (rather than simply omitted) so it's visibly a choice.
    breakdown.surface = 0.0;

    return breakdown;
}

void Ewald::computeForces(const core::ParticleData& particles, const core::Lattice& lattice,
                           const core::NeighborList& neighbors, Forces& forcesOut) const {
    validateChargeNeutrality(particles);
    const ExclusionSet excluded(exclusions_);

    // --- Real-space force: F_i = (dV/dr / r) * d = -factor * d, same sign
    // convention as LennardJones::computeForces. ---
    for (const auto& pair : neighbors.pairs()) {
        if (excluded.contains(pair.i, pair.j)) {
            continue;
        }
        const Vec3 d = lattice.minimumImageDisplacement(position(particles, pair.i),
                                                          position(particles, pair.j));
        const double r = norm(d);
        if (r > realSpaceCutoff_) {
            continue;
        }
        const double factor = kCoulombConstant * particles.charge[pair.i] *
                               particles.charge[pair.j] * realSpaceForceOverR(r, alpha_);
        forcesOut.fx[pair.i] -= factor * d[0];
        forcesOut.fy[pair.i] -= factor * d[1];
        forcesOut.fz[pair.i] -= factor * d[2];
        forcesOut.fx[pair.j] += factor * d[0];
        forcesOut.fy[pair.j] += factor * d[1];
        forcesOut.fz[pair.j] += factor * d[2];
    }

    // --- Exclusion force correction: negative gradient of the subtracted
    // -q_i*q_j*erf(alpha*r)/r energy term. Same F_i = -factor*d
    // convention. ---
    for (const auto& [i, j] : exclusions_) {
        const Vec3 d =
            lattice.minimumImageDisplacement(position(particles, i), position(particles, j));
        const double r = norm(d);
        const double factor =
            -kCoulombConstant * particles.charge[i] * particles.charge[j] * exclusionForceOverR(r, alpha_);
        forcesOut.fx[i] -= factor * d[0];
        forcesOut.fy[i] -= factor * d[1];
        forcesOut.fz[i] -= factor * d[2];
        forcesOut.fx[j] += factor * d[0];
        forcesOut.fy[j] += factor * d[1];
        forcesOut.fz[j] += factor * d[2];
    }

    // --- Reciprocal-space force:
    // F_i = -(4*pi/V) * kCoulomb * q_i * sum_k weight(k) *
    //       [Im(k)*cos(k.r_i) - Re(k)*sin(k.r_i)] * k
    // (negative gradient of U_reciprocal w.r.t. r_i; weight(k) is the same
    // exp(-k^2/4a^2)/k^2 factor used in the energy, and Re(k)/Im(k) are
    // the same structure-factor components — computed once per k and
    // reused, so energy and force can never disagree about what S(k) is).
    // ---
    const double volume = std::abs(lattice.volume());
    const auto reciprocalTerms = buildReciprocalTerms(lattice, alpha_, kMax_);
    const double prefactor = kCoulombConstant * 4.0 * std::numbers::pi / volume;
    for (const auto& term : reciprocalTerms) {
        const auto s = structureFactor(particles, term.k);
        for (std::size_t i = 0; i < particles.size(); ++i) {
            const double kr = dot(term.k, position(particles, i));
            const double coeff =
                term.weight * (s.im * std::cos(kr) - s.re * std::sin(kr)) * particles.charge[i];
            const double scaleFactor = -prefactor * coeff;
            forcesOut.fx[i] += scaleFactor * term.k[0];
            forcesOut.fy[i] += scaleFactor * term.k[1];
            forcesOut.fz[i] += scaleFactor * term.k[2];
        }
    }
}

} // namespace aleator::forcefield
