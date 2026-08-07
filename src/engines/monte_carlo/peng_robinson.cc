#include "engines/monte_carlo/peng_robinson.hpp"

#include <cmath>
#include <numbers>
#include <stdexcept>
#include <string>

namespace aleator::engines {

PengRobinson::PengRobinson(PengRobinsonSubstance substance) : substance_(substance) {
    if (!(substance_.criticalTemperatureKelvin > 0.0)) {
        throw std::invalid_argument("PengRobinson: criticalTemperatureKelvin must be positive");
    }
    if (!(substance_.criticalPressurePascal > 0.0)) {
        throw std::invalid_argument("PengRobinson: criticalPressurePascal must be positive");
    }
    const double tc = substance_.criticalTemperatureKelvin;
    const double pc = substance_.criticalPressurePascal;
    const double r = kGasConstant;
    a_ = 0.45724 * r * r * tc * tc / pc;
    b_ = 0.07780 * r * tc / pc;
    const double omega = substance_.acentricFactor;
    kappa_ = 0.37464 + 1.54226 * omega - 0.26992 * omega * omega;
}

namespace {
double alphaOfT(double temperatureKelvin, double criticalTemperatureKelvin, double kappa) {
    const double term = 1.0 + kappa * (1.0 - std::sqrt(temperatureKelvin / criticalTemperatureKelvin));
    return term * term;
}
} // namespace

double PengRobinson::pressurePascal(double temperatureKelvin, double molarVolumeM3PerMol) const {
    const double v = molarVolumeM3PerMol;
    const double aAlpha = a_ * alphaOfT(temperatureKelvin, substance_.criticalTemperatureKelvin, kappa_);
    return (kGasConstant * temperatureKelvin) / (v - b_) -
           aAlpha / (v * v + 2.0 * b_ * v - b_ * b_);
}

double PengRobinson::compressibilityFactor(double temperatureKelvin, double pressurePascal) const {
    if (!(temperatureKelvin > 0.0)) {
        throw std::invalid_argument("PengRobinson: temperatureKelvin must be positive");
    }
    if (!(pressurePascal > 0.0)) {
        throw std::invalid_argument("PengRobinson: pressurePascal must be positive");
    }
    const double aAlpha =
        a_ * alphaOfT(temperatureKelvin, substance_.criticalTemperatureKelvin, kappa_);
    const double r = kGasConstant;
    const double bigA = aAlpha * pressurePascal / (r * r * temperatureKelvin * temperatureKelvin);
    const double bigB = b_ * pressurePascal / (r * temperatureKelvin);

    // Z^3 - (1-B)Z^2 + (A - 3B^2 - 2B)Z - (AB - B^2 - B^3) = 0
    const double c2 = -(1.0 - bigB);
    const double c1 = bigA - 3.0 * bigB * bigB - 2.0 * bigB;
    const double c0 = -(bigA * bigB - bigB * bigB - bigB * bigB * bigB);

    auto f = [&](double z) { return z * z * z + c2 * z * z + c1 * z + c0; };
    auto fPrime = [&](double z) { return 3.0 * z * z + 2.0 * c2 * z + c1; };

    // Newton-Raphson from the ideal-gas root Z=1 -- the vapor/gas-phase
    // root, correct for the T > Tc (supercritical) states this codebase's
    // GCMC use targets. Not a general three-real-root cubic classifier
    // (that would matter for subcritical vapor-liquid equilibrium, out of
    // scope here — see class doc comment).
    double z = 1.0;
    bool converged = false;
    for (int iter = 0; iter < 100; ++iter) {
        const double fz = f(z);
        if (std::abs(fz) < 1e-14) {
            converged = true;
            break;
        }
        const double fpz = fPrime(z);
        if (fpz == 0.0) {
            break;
        }
        z -= fz / fpz;
        if (!(z > bigB)) {
            // Stepped into the unphysical Z<=B region (ln(Z-B) undefined);
            // back off toward the ideal-gas side rather than diverge.
            z = bigB + 0.5 * (1.0 - bigB > 0 ? 1.0 - bigB : 1.0);
        }
    }
    if (!converged && std::abs(f(z)) > 1e-8) {
        throw std::runtime_error(
            "PengRobinson::compressibilityFactor: Newton-Raphson did not converge for T=" +
            std::to_string(temperatureKelvin) + " K, P=" + std::to_string(pressurePascal) + " Pa");
    }
    return z;
}

double PengRobinson::fugacityCoefficient(double temperatureKelvin, double pressurePascal) const {
    const double z = compressibilityFactor(temperatureKelvin, pressurePascal);
    const double aAlpha =
        a_ * alphaOfT(temperatureKelvin, substance_.criticalTemperatureKelvin, kappa_);
    const double r = kGasConstant;
    const double bigA = aAlpha * pressurePascal / (r * r * temperatureKelvin * temperatureKelvin);
    const double bigB = b_ * pressurePascal / (r * temperatureKelvin);
    const double sqrt2 = std::numbers::sqrt2;

    const double lnPhi = (z - 1.0) - std::log(z - bigB) -
                          (bigA / (2.0 * sqrt2 * bigB)) *
                              std::log((z + (1.0 + sqrt2) * bigB) / (z + (1.0 - sqrt2) * bigB));
    return std::exp(lnPhi);
}

double PengRobinson::fugacityPascal(double temperatureKelvin, double pressurePascal) const {
    return fugacityCoefficient(temperatureKelvin, pressurePascal) * pressurePascal;
}

PengRobinsonMixture::PengRobinsonMixture(std::vector<PengRobinsonSubstance> substances,
                                          std::vector<std::vector<double>> binaryInteractionKij)
    : substances_(std::move(substances)) {
    if (substances_.size() < 2) {
        throw std::invalid_argument(
            "PengRobinsonMixture: needs at least 2 substances (use PengRobinson for a pure "
            "component)");
    }
    const std::size_t n = substances_.size();
    a_.reserve(n);
    b_.reserve(n);
    kappa_.reserve(n);
    for (const auto& substance : substances_) {
        if (!(substance.criticalTemperatureKelvin > 0.0)) {
            throw std::invalid_argument(
                "PengRobinsonMixture: criticalTemperatureKelvin must be positive");
        }
        if (!(substance.criticalPressurePascal > 0.0)) {
            throw std::invalid_argument(
                "PengRobinsonMixture: criticalPressurePascal must be positive");
        }
        const double tc = substance.criticalTemperatureKelvin;
        const double pc = substance.criticalPressurePascal;
        const double r = PengRobinson::kGasConstant;
        a_.push_back(0.45724 * r * r * tc * tc / pc);
        b_.push_back(0.07780 * r * tc / pc);
        const double omega = substance.acentricFactor;
        kappa_.push_back(0.37464 + 1.54226 * omega - 0.26992 * omega * omega);
    }

    if (binaryInteractionKij.empty()) {
        kij_.assign(n, std::vector<double>(n, 0.0));
    } else {
        if (binaryInteractionKij.size() != n) {
            throw std::invalid_argument(
                "PengRobinsonMixture: binaryInteractionKij must be substances.size() x "
                "substances.size()");
        }
        for (std::size_t i = 0; i < n; ++i) {
            if (binaryInteractionKij[i].size() != n) {
                throw std::invalid_argument(
                    "PengRobinsonMixture: binaryInteractionKij must be substances.size() x "
                    "substances.size()");
            }
            if (binaryInteractionKij[i][i] != 0.0) {
                throw std::invalid_argument(
                    "PengRobinsonMixture: binaryInteractionKij diagonal must be zero (k_ii = 0)");
            }
            for (std::size_t j = i + 1; j < n; ++j) {
                if (binaryInteractionKij[i][j] != binaryInteractionKij[j][i]) {
                    throw std::invalid_argument(
                        "PengRobinsonMixture: binaryInteractionKij must be symmetric (k_ij = k_ji)");
                }
            }
        }
        kij_ = std::move(binaryInteractionKij);
    }
}

namespace {

void validateMoleFractions(const std::vector<double>& moleFractions, std::size_t expectedSize) {
    if (moleFractions.size() != expectedSize) {
        throw std::invalid_argument(
            "PengRobinsonMixture: moleFractions.size() must equal the number of substances");
    }
    double sum = 0.0;
    for (double x : moleFractions) {
        if (!(x >= 0.0)) {
            throw std::invalid_argument("PengRobinsonMixture: every mole fraction must be >= 0");
        }
        sum += x;
    }
    constexpr double kSumTolerance = 1e-9;
    if (std::abs(sum - 1.0) > kSumTolerance) {
        throw std::invalid_argument("PengRobinsonMixture: mole fractions must sum to 1 (got " +
                                     std::to_string(sum) + ")");
    }
}

} // namespace

double PengRobinsonMixture::compressibilityFactor(double temperatureKelvin,
                                                    double totalPressurePascal,
                                                    const std::vector<double>& moleFractions) const {
    if (!(temperatureKelvin > 0.0)) {
        throw std::invalid_argument("PengRobinsonMixture: temperatureKelvin must be positive");
    }
    if (!(totalPressurePascal > 0.0)) {
        throw std::invalid_argument("PengRobinsonMixture: totalPressurePascal must be positive");
    }
    validateMoleFractions(moleFractions, substances_.size());
    const std::size_t n = substances_.size();

    // a_i(T) = a_i * alpha_i(T), each component's own temperature-dependent
    // pure-substance PR "a" -- see class doc comment.
    std::vector<double> aAlpha(n);
    for (std::size_t i = 0; i < n; ++i) {
        aAlpha[i] = a_[i] * alphaOfT(temperatureKelvin, substances_[i].criticalTemperatureKelvin,
                                       kappa_[i]);
    }

    double aMix = 0.0;
    double bMix = 0.0;
    for (std::size_t i = 0; i < n; ++i) {
        bMix += moleFractions[i] * b_[i];
        for (std::size_t j = 0; j < n; ++j) {
            const double aij = std::sqrt(aAlpha[i] * aAlpha[j]) * (1.0 - kij_[i][j]);
            aMix += moleFractions[i] * moleFractions[j] * aij;
        }
    }

    const double r = PengRobinson::kGasConstant;
    const double bigA =
        aMix * totalPressurePascal / (r * r * temperatureKelvin * temperatureKelvin);
    const double bigB = bMix * totalPressurePascal / (r * temperatureKelvin);

    // Identical cubic form (and identical Newton-Raphson-from-Z=1 solve,
    // gas/vapor-phase root only) as PengRobinson::compressibilityFactor,
    // parameterized by the MIXTURE's A/B rather than a pure substance's.
    const double c2 = -(1.0 - bigB);
    const double c1 = bigA - 3.0 * bigB * bigB - 2.0 * bigB;
    const double c0 = -(bigA * bigB - bigB * bigB - bigB * bigB * bigB);

    auto f = [&](double z) { return z * z * z + c2 * z * z + c1 * z + c0; };
    auto fPrime = [&](double z) { return 3.0 * z * z + 2.0 * c2 * z + c1; };

    double z = 1.0;
    bool converged = false;
    for (int iter = 0; iter < 100; ++iter) {
        const double fz = f(z);
        if (std::abs(fz) < 1e-14) {
            converged = true;
            break;
        }
        const double fpz = fPrime(z);
        if (fpz == 0.0) {
            break;
        }
        z -= fz / fpz;
        if (!(z > bigB)) {
            z = bigB + 0.5 * (1.0 - bigB > 0 ? 1.0 - bigB : 1.0);
        }
    }
    if (!converged && std::abs(f(z)) > 1e-8) {
        throw std::runtime_error(
            "PengRobinsonMixture::compressibilityFactor: Newton-Raphson did not converge for T=" +
            std::to_string(temperatureKelvin) + " K, P=" + std::to_string(totalPressurePascal) +
            " Pa");
    }
    return z;
}

std::vector<double> PengRobinsonMixture::fugacitiesPascal(
    double temperatureKelvin, double totalPressurePascal,
    const std::vector<double>& moleFractions) const {
    validateMoleFractions(moleFractions, substances_.size());
    const std::size_t n = substances_.size();
    const double z = compressibilityFactor(temperatureKelvin, totalPressurePascal, moleFractions);

    std::vector<double> aAlpha(n);
    for (std::size_t i = 0; i < n; ++i) {
        aAlpha[i] = a_[i] * alphaOfT(temperatureKelvin, substances_[i].criticalTemperatureKelvin,
                                       kappa_[i]);
    }
    double aMix = 0.0;
    double bMix = 0.0;
    for (std::size_t i = 0; i < n; ++i) {
        bMix += moleFractions[i] * b_[i];
        for (std::size_t j = 0; j < n; ++j) {
            aMix += moleFractions[i] * moleFractions[j] * std::sqrt(aAlpha[i] * aAlpha[j]) *
                    (1.0 - kij_[i][j]);
        }
    }

    const double r = PengRobinson::kGasConstant;
    const double bigA =
        aMix * totalPressurePascal / (r * r * temperatureKelvin * temperatureKelvin);
    const double bigB = bMix * totalPressurePascal / (r * temperatureKelvin);
    const double sqrt2 = std::numbers::sqrt2;

    std::vector<double> fugacities(n);
    for (std::size_t i = 0; i < n; ++i) {
        double sumJ = 0.0;
        for (std::size_t j = 0; j < n; ++j) {
            sumJ += moleFractions[j] * std::sqrt(aAlpha[i] * aAlpha[j]) * (1.0 - kij_[i][j]);
        }
        const double lnPhiI =
            (b_[i] / bMix) * (z - 1.0) - std::log(z - bigB) -
            (bigA / (2.0 * sqrt2 * bigB)) * ((2.0 / aMix) * sumJ - b_[i] / bMix) *
                std::log((z + (1.0 + sqrt2) * bigB) / (z + (1.0 - sqrt2) * bigB));
        const double phiI = std::exp(lnPhiI);
        fugacities[i] = moleFractions[i] * phiI * totalPressurePascal;
    }
    return fugacities;
}

} // namespace aleator::engines
