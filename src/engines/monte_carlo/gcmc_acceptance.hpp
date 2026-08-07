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

/// Unclamped ratio for an identity-swap move: converting one existing
/// molecule of species A (countABefore molecules of A present before the
/// move) into a molecule of species B (countBBefore molecules of B present
/// before the move), in place.
///
/// Derivation (this is "where the classic errors live" per the mixture
/// milestone that added this function — worth writing out in full rather
/// than asserting the result). The multi-species grand-canonical weight is
///   pi(N_A, N_B, ...) ~ prod_i [z_i^{N_i} V^{N_i} / (N_i! Lambda_i^{3N_i})] * exp(-beta*U)
/// with each species' activity z_i = beta * Lambda_i^3 * f_i (see
/// insertionRatio's doc comment for why Lambda_i cancels completely within
/// species i's own factor: z_i's Lambda_i^3 exactly cancels the
/// Lambda_i^{-3N_i} in that species' own normalization, regardless of what
/// any OTHER species is doing — there is no cross term between species in an
/// ideal, non-interacting reservoir). For the transition (N_A, N_B) ->
/// (N_A-1, N_B+1):
///   pi(N_A-1,N_B+1) / pi(N_A,N_B)
///     = (z_B/z_A) * (N_A! / (N_A-1)!) * (N_B! / (N_B+1)!)
///     = (z_B/z_A) * N_A / (N_B+1)
/// and z_B/z_A = f_B/f_A exactly (each species' own Lambda_i cancels within
/// its own z_i, per above — there is no leftover (Lambda_A/Lambda_B)^3 mass
/// factor, even though A and B are physically different molecules with
/// different masses and therefore different Lambda). So:
///   R_swap(A->B) = (f_B/f_A) * (countABefore / (countBBefore+1)) * exp(-deltaU/T)
/// with NO volume or temperature prefactor beyond what's inside deltaU's
/// exp() — V/T cancel exactly, confirmed independently by composing this as
/// "delete A then insert B at the same position": deletionRatio(fA, V,
/// countABefore, T, dU_A) * insertionRatio(fB, V, countBBefore, T, dU_B) =
/// [countABefore/(fA*V/T)] * [(fB*V/T)/(countBBefore+1)] =
/// (fB*countABefore)/(fA*(countBBefore+1)), the V/T cancelling identically —
/// the same answer by a second, independent route.
///
/// Reciprocal-identity check (tests/validation/test_gcmc_detailed_balance.cc):
/// swapRatio(fA,fB,countABefore,countBBefore,T,dU) *
/// swapRatio(fB,fA,countBBefore+1,countABefore-1,T,-dU) == 1 exactly, on the
/// reverse transition (N_A-1,N_B+1) -> (N_A,N_B).
[[nodiscard]] inline double swapRatio(double fugacityInternalFrom, double fugacityInternalTo,
                                       std::size_t countFromBefore, std::size_t countToBefore,
                                       double temperatureKelvin, double deltaUKelvin) {
    return (fugacityInternalTo / fugacityInternalFrom) *
           (static_cast<double>(countFromBefore) / static_cast<double>(countToBefore + 1)) *
           std::exp(-deltaUKelvin / temperatureKelvin);
}

/// min(1, ratio) — the acceptance PROBABILITY actually used by the engine;
/// the *Ratio functions above return the unclamped value because the
/// detailed-balance identity (product of forward/reverse ratios == 1) only
/// holds before clamping.
[[nodiscard]] inline double clampToProbability(double ratio) { return ratio < 1.0 ? ratio : 1.0; }

} // namespace aleator::engines::gcmc
