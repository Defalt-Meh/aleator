#pragma once

#include <cmath>
#include <cstddef>

namespace aleator::engines::gcmc {

/// Unclamped Metropolis-style acceptance RATIO for inserting one molecule
/// into a system that currently has `countBefore` molecules of this
/// species. Derived from the grand-canonical configurational weight
/// pi(N) ~ (z^N V^N)/(N! Lambda^{3N}) * exp(-beta*U(N)), activity z =
/// beta*Lambda^3*f (f = fugacity) — under which the thermal de Broglie
/// wavelength Lambda cancels completely, so it never needs to be computed
/// at all:
///   R_insert = (f*V/T) / (countBefore + 1) * exp(-deltaU/T)
/// This codebase's internal energy unit is already K = energy/k_B, so
/// beta*E is simply E[K]/T[K] — no separate Boltzmann constant appears.
/// `fugacityInternal` must be in this codebase's internal pressure/energy-
/// density convention, K/Angstrom^3 (see MonteCarloEngine's doc comment
/// for the Pa -> K/Angstrom^3 conversion), `volume` in Angstrom^3.
///
/// See tests/validation/test_gcmc_detailed_balance.cc: insertionRatio and
/// deletionRatio are exact reciprocals of each other on the same pair of
/// states, which — under the standard min(1,R)/min(1,1/R) Metropolis
/// acceptance prescription — is detailed balance for this move pair,
/// checked there as a runtime-verifiable algebraic identity rather than a
/// code-review claim.
[[nodiscard]] inline double insertionRatio(double fugacityInternal, double volume,
                                            std::size_t countBefore, double temperatureKelvin,
                                            double deltaUKelvin) {
    return (fugacityInternal * volume / temperatureKelvin) /
           static_cast<double>(countBefore + 1) * std::exp(-deltaUKelvin / temperatureKelvin);
}

/// Unclamped ratio for deleting one (uniformly randomly chosen) molecule
/// from a system that currently has `countBefore` molecules of this
/// species:
///   R_delete = countBefore / (f*V/T) * exp(-deltaU/T)
[[nodiscard]] inline double deletionRatio(double fugacityInternal, double volume,
                                           std::size_t countBefore, double temperatureKelvin,
                                           double deltaUKelvin) {
    return static_cast<double>(countBefore) / (fugacityInternal * volume / temperatureKelvin) *
           std::exp(-deltaUKelvin / temperatureKelvin);
}

/// Unclamped ratio for a symmetric-proposal move (translation, rotation)
/// that doesn't change the particle count: plain Metropolis,
/// R = exp(-deltaU/T).
[[nodiscard]] inline double metropolisRatio(double temperatureKelvin, double deltaUKelvin) {
    return std::exp(-deltaUKelvin / temperatureKelvin);
}

/// min(1, ratio) — the acceptance PROBABILITY actually used by the engine;
/// the *Ratio functions above return the unclamped value because the
/// detailed-balance identity (product of forward/reverse ratios == 1) only
/// holds before clamping.
[[nodiscard]] inline double clampToProbability(double ratio) { return ratio < 1.0 ? ratio : 1.0; }

} // namespace aleator::engines::gcmc
