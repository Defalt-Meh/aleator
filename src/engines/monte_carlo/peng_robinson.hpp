#pragma once

#include <vector>

namespace aleator::engines {

/// Pure-substance critical properties needed by the Peng-Robinson EOS.
/// Methane's values (CODATA/NIST-consistent, cross-checked against the
/// original Peng & Robinson (1976) parameterization): critical temperature
/// 190.6 K, critical pressure 4.599 MPa, acentric factor 0.011.
struct PengRobinsonSubstance {
    double criticalTemperatureKelvin;
    double criticalPressurePascal;
    double acentricFactor;
};

inline constexpr PengRobinsonSubstance kMethane{190.6, 4.599e6, 0.011};

/// Peng-Robinson (1976) cubic equation of state for a pure substance:
///   P = RT/(V-b) - a*alpha(T) / (V^2 + 2bV - b^2)
/// used here for one purpose: converting an external gas-phase pressure
/// into the fugacity that appears in the GCMC insertion/deletion acceptance
/// criteria (engines/monte_carlo/gcmc_acceptance.hpp) — methane at the
/// pressures/temperatures relevant to adsorption isotherms deviates
/// non-negligibly from ideal-gas behavior, and using fugacity rather than
/// raw pressure there is what keeps the chemical-potential term physically
/// correct.
///
/// a = 0.45724 * R^2 * Tc^2 / Pc, b = 0.07780 * R * Tc / Pc, kappa =
/// 0.37464 + 1.54226*omega - 0.26992*omega^2, alpha(T) = (1 + kappa*(1 -
/// sqrt(T/Tc)))^2 — the coefficients 0.45724/0.07780/0.37464/1.54226/
/// 0.26992 are Peng & Robinson's own published constants (Peng, D.Y.;
/// Robinson, D.B. "A New Two-Constant Equation of State." Ind. Eng. Chem.
/// Fundam. 1976, 15 (1), 59-64), not re-derived or invented here.
///
/// All quantities SI: T in K, P in Pa, molar volume in m^3/mol. Internal
/// only to this class — the GCMC engine converts the resulting fugacity
/// (Pa) to this codebase's internal K/Angstrom^3 energy-density convention
/// at the point of use.
class PengRobinson {
public:
    explicit PengRobinson(PengRobinsonSubstance substance);

    /// Molar gas constant, R = N_A * k_B (SI 2019 exact via the exact
    /// Avogadro and Boltzmann constants): 8.31446261815324 J/(mol*K).
    static constexpr double kGasConstant = 8.31446261815324;

    /// Compressibility factor Z = PV/(RT) for the vapor/gas-phase root of
    /// the PR cubic, found by Newton-Raphson from the ideal-gas starting
    /// point Z=1. Correct for the supercritical/gas states this codebase's
    /// GCMC use targets (methane at typical adsorption-isotherm conditions,
    /// T > Tc); throws std::runtime_error if Newton-Raphson fails to
    /// converge rather than returning a plausible-looking wrong root.
    [[nodiscard]] double compressibilityFactor(double temperatureKelvin,
                                                double pressurePascal) const;

    /// Fugacity coefficient phi = f/P (dimensionless), from the standard
    /// closed-form PR result:
    ///   ln(phi) = Z - 1 - ln(Z - B) -
    ///             (A / (2*sqrt(2)*B)) * ln[(Z + (1+sqrt2)B) / (Z + (1-sqrt2)B)]
    /// with A = a*alpha(T)*P/(R^2 T^2), B = b*P/(R T).
    [[nodiscard]] double fugacityCoefficient(double temperatureKelvin, double pressurePascal) const;

    /// f = phi * P, in Pa.
    [[nodiscard]] double fugacityPascal(double temperatureKelvin, double pressurePascal) const;

    /// P(T, V) directly from the EOS (not through the cubic solver) —
    /// exposed for validating that a/b/alpha(T) genuinely reproduce a
    /// critical point at (Tc, Pc): dP/dV and d^2P/dV^2 must both vanish
    /// there (tests/validation/test_peng_robinson.cc checks this via finite
    /// difference, independent of compressibilityFactor()'s Newton solve).
    [[nodiscard]] double pressurePascal(double temperatureKelvin, double molarVolumeM3PerMol) const;

    /// The EOS covolume parameter b (m^3/mol), exposed for testing: at
    /// T=Tc, P=Pc the dimensionless B = b*Pc/(R*Tc) is EXACTLY 0.07780 by
    /// construction (b is defined as 0.07780*R*Tc/Pc), which combined with
    /// alpha(Tc)=1 makes the critical Z a closed-form root of the cubic
    /// (Z - Zc)^3 = 0 with Zc = (1-B_c)/3 — letting the critical-point
    /// self-consistency test check dP/dV=d^2P/dV^2=0 there directly,
    /// without going through compressibilityFactor()'s Newton solve (which
    /// is numerically ill-conditioned exactly at a triple root, since f'
    /// vanishes there too).
    [[nodiscard]] double bParameter() const noexcept { return b_; }

private:
    PengRobinsonSubstance substance_;
    double a_;
    double b_;
    double kappa_;
};

/// Peng-Robinson (1976) cubic EOS for a MIXTURE of substances, via the
/// standard van der Waals one-fluid mixing rules (Peng & Robinson's own
/// paper already extends the pure-substance EOS this way; the specific
/// closed-form per-component fugacity-coefficient expression used below
/// follows the presentation in Reid/Prausnitz/Poling, "The Properties of
/// Gases and Liquids," and is standard across the mixture-GCMC/screening
/// literature, not derived from scratch here):
///
///   a_ij(T) = sqrt(a_i(T) * a_j(T)) * (1 - k_ij)     [a_i(T) = a_i * alpha_i(T),
///                                                      the pure component's own
///                                                      temperature-dependent PR "a"]
///   a_mix(T) = sum_i sum_j x_i x_j a_ij(T)
///   b_mix = sum_i x_i b_i
///
/// with k_ij the (symmetric, zero-diagonal) binary interaction parameter
/// between components i and j. THIS CLASS DEFAULTS k_ij = 0 for every pair
/// UNLESS EXPLICITLY GIVEN -- an honest, documented simplification (no
/// regressed binary interaction data is bundled with this codebase), not a
/// silent assumption: k_ij = 0 is itself a real, named approximation
/// (equivalent to assuming the unlike-pair interaction is exactly the
/// geometric mean of the like-pair interactions), and callers with real
/// k_ij data for their system should supply it via the constructor.
///
/// Per-component fugacity coefficient in the mixture, at the mixture's own
/// Z (found from A_mix/B_mix via the identical cubic PengRobinson::
/// compressibilityFactor solves for a pure substance):
///   ln(phi_i) = (b_i/b_mix)*(Z-1) - ln(Z-B_mix)
///     - (A_mix/(2*sqrt2*B_mix)) * [(2/a_mix)*sum_j(x_j*a_ij) - b_i/b_mix]
///       * ln[(Z+(1+sqrt2)*B_mix)/(Z+(1-sqrt2)*B_mix)]
/// and f_i = x_i * phi_i * P (the component's fugacity AT its mole fraction
/// in the mixture, i.e. its partial fugacity -- this is what
/// MonteCarloEngine's per-species GCMC insertion/deletion/swap acceptance
/// (engines/monte_carlo/gcmc_acceptance.hpp) needs for species i, not the
/// pure-component fugacity PengRobinson::fugacityPascal would give at the
/// same total pressure).
class PengRobinsonMixture {
public:
    /// `binaryInteractionKij`, if given, must be substances.size() x
    /// substances.size(), symmetric, zero diagonal -- checked at
    /// construction. Empty (the default) means every k_ij = 0 (see class
    /// doc comment).
    explicit PengRobinsonMixture(std::vector<PengRobinsonSubstance> substances,
                                  std::vector<std::vector<double>> binaryInteractionKij = {});

    /// Mixture compressibility factor Z at the given T, total P, and mole
    /// fractions (must sum to 1, checked). Vapor/gas-phase root only, same
    /// caveats as PengRobinson::compressibilityFactor.
    [[nodiscard]] double compressibilityFactor(double temperatureKelvin,
                                                double totalPressurePascal,
                                                const std::vector<double>& moleFractions) const;

    /// Per-component fugacity (Pa), f_i = x_i * phi_i * P, at the given T,
    /// total P, and mole fractions.
    [[nodiscard]] std::vector<double> fugacitiesPascal(
        double temperatureKelvin, double totalPressurePascal,
        const std::vector<double>& moleFractions) const;

private:
    std::vector<PengRobinsonSubstance> substances_;
    std::vector<double> a_;     // per-component pure PR "a" (same formula as PengRobinson)
    std::vector<double> b_;     // per-component pure PR "b"
    std::vector<double> kappa_; // per-component PR kappa
    std::vector<std::vector<double>> kij_;
};

} // namespace aleator::engines
