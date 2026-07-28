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

} // namespace aleator::engines
