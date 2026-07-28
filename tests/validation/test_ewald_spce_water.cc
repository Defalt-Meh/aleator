// Validation anchor (CLAUDE.md milestone: "NIST SRSW SPC/E water Ewald
// reference energies"). Uses real published NIST SRSW SPC/E water sample
// configurations (see data/nist_spce/PROVENANCE.md; not invented) and
// compares each individual Ewald term (real, reciprocal/"fourier", self,
// intramolecular exclusion/"intra") against NIST's published per-term
// breakdown -- not just the total -- since E_self and E_intra are each
// O(1e6-1e7) K and nearly cancel, so a total-only comparison would hide a
// term-coefficient bug behind that cancellation. This is exactly the
// failure mode CLAUDE.md #10 warns about ("Ewald missing the self-
// interaction or exclusion correction -- smooth, plausible, wrong").

#include <catch2/catch_test_macros.hpp>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdio>
#include <filesystem>
#include <fstream>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include "core/geometry/lattice.hpp"
#include "core/math/particle_data.hpp"
#include "core/neighbor/verlet_list.hpp"
#include "forcefield/electrostatics/ewald.hpp"

namespace {

using aleator::core::Lattice;
using aleator::core::ParticleData;
using aleator::core::VerletList;
using aleator::forcefield::Ewald;
using aleator::forcefield::EwaldEnergyBreakdown;

// SPC/E model parameters, verbatim from the NIST SRSW page (see
// data/nist_spce/PROVENANCE.md): H = +q, O = -2q.
constexpr double kSpceCharge = 0.42380; // e

struct SpceConfig {
    Lattice lattice;
    ParticleData particles;
    std::vector<std::pair<std::size_t, std::size_t>> exclusions;
};

/// Parses a NIST SRSW SPC/E sample configuration file (see
/// data/nist_spce/metadata.README): line 1 = box x/y/z extent, line 2 =
/// molecule count M, then 3*M rows of "id x y z type". Every bundled file
/// groups atoms O,H,H per molecule in that fixed order (verified by hand
/// for all four -- see the shell check referenced in this milestone's
/// notes); this parser checks that assumption explicitly per atom rather
/// than trusting it silently, throwing a clear error if it ever doesn't
/// hold, per CLAUDE.md #7 ("no stub/silent-wrong-parse").
SpceConfig loadSpceConfig(const std::filesystem::path& path) {
    std::ifstream in(path);
    if (!in) {
        throw std::runtime_error("could not open NIST SPC/E config file: " + path.string());
    }

    double lx = 0.0;
    double ly = 0.0;
    double lz = 0.0;
    in >> lx >> ly >> lz;

    std::size_t moleculeCount = 0;
    in >> moleculeCount;
    const std::size_t atomCount = 3 * moleculeCount;

    ParticleData particles;
    particles.resize(atomCount);
    std::vector<std::pair<std::size_t, std::size_t>> exclusions;
    exclusions.reserve(3 * moleculeCount);

    for (std::size_t k = 0; k < atomCount; ++k) {
        std::size_t id = 0;
        double x = 0.0;
        double y = 0.0;
        double z = 0.0;
        std::string type;
        in >> id >> x >> y >> z >> type;
        if (!in) {
            throw std::runtime_error("malformed NIST SPC/E config file (short read) at atom " +
                                      std::to_string(k) + ": " + path.string());
        }
        const bool expectOxygen = (k % 3 == 0);
        if (expectOxygen != (type == "O")) {
            throw std::runtime_error(
                "NIST SPC/E config file does not follow the assumed O,H,H per-molecule "
                "ordering at atom " +
                std::to_string(k) + ": " + path.string());
        }
        particles.x[k] = x;
        particles.y[k] = y;
        particles.z[k] = z;
        particles.species[k] = expectOxygen ? 0 : 1;
        particles.charge[k] = expectOxygen ? -2.0 * kSpceCharge : kSpceCharge;
    }

    // All 3 pairwise combinations within each rigid 3-site molecule are
    // excluded from direct Coulomb (accounted for instead by U_exclusion) --
    // matches NIST's "E_intra" sum over all site pairs within a molecule.
    for (std::size_t m = 0; m < moleculeCount; ++m) {
        const std::size_t o = 3 * m;
        const std::size_t h1 = 3 * m + 1;
        const std::size_t h2 = 3 * m + 2;
        exclusions.emplace_back(o, h1);
        exclusions.emplace_back(o, h2);
        exclusions.emplace_back(h1, h2);
    }

    std::array<std::array<double, 3>, 3> m{};
    m[0] = {lx, 0.0, 0.0};
    m[1] = {0.0, ly, 0.0};
    m[2] = {0.0, 0.0, lz};
    return {Lattice(m), std::move(particles), std::move(exclusions)};
}

std::filesystem::path spceDataFile(const std::string& name) {
    return std::filesystem::path(ALEATOR_TEST_DATA_DIR) / "data" / "nist_spce" / name;
}

struct SpceReference {
    std::string configFile;
    std::size_t moleculeCount;
    double eReal;
    double eFourier;
    double eSelf;
    double eIntra;
};

// Transcribed from the NIST SRSW SPC/E water reference-calculations table
// (10 Angstrom real-space cutoff); see data/nist_spce/PROVENANCE.md for the
// exact source URL. Energies are per k_B, in K.
const std::vector<SpceReference>& nistSpceReferenceData() {
    static const std::vector<SpceReference> data{
        {"spce_sample_config_periodic1.txt", 100, -5.58889e+05, 6.27009e+03, -2.84469e+06,
         2.80999e+06},
        {"spce_sample_config_periodic2.txt", 200, -1.19295e+06, 6.03495e+03, -5.68938e+06,
         5.61998e+06},
        {"spce_sample_config_periodic3.txt", 300, -1.96297e+06, 5.24461e+03, -8.53407e+06,
         8.42998e+06},
        {"spce_sample_config_periodic4.txt", 750, -3.57226e+06, 7.58785e+03, -1.42235e+07,
         1.41483e+07},
    };
    return data;
}

// NIST's published values carry 6 significant figures and were computed
// with CODATA 2010 physical constants (~1.5e-6 relative different from
// this codebase's CODATA 2018/SI-2019 Ewald::kCoulombConstant) -- see
// data/nist_spce/PROVENANCE.md. A flat 1e-8 (CLAUDE.md's LJ-energy figure)
// is not achievable against a 6-sig-fig reference; this mirrors the
// LJ-NIST test's tolerance calibration (test_lennard_jones_nist.cc), not a
// loosened-to-pass number (CLAUDE.md invariant #1).
constexpr double kRelativeTolerance = 5e-4;

} // namespace

TEST_CASE("Ewald matches the published NIST SRSW SPC/E water reference energy breakdown",
          "[validation][electrostatics]") {
    for (const auto& ref : nistSpceReferenceData()) {
        const SpceConfig config = loadSpceConfig(spceDataFile(ref.configFile));
        REQUIRE(config.particles.size() == 3 * ref.moleculeCount);

        const auto& latticeMatrix = config.lattice.matrix();
        const double minBoxLength = std::min(
            {latticeMatrix[0][0], latticeMatrix[1][1], latticeMatrix[2][2]});
        const double alpha = 5.6 / minBoxLength; // NIST's exact stated convention
        const double realSpaceCutoff = 10.0;
        const int kMax = 5;

        VerletList neighbors(/*skinDistance=*/0.1);
        neighbors.build(config.particles, config.lattice, realSpaceCutoff);

        const Ewald ewald(alpha, realSpaceCutoff, kMax, config.exclusions);
        const EwaldEnergyBreakdown breakdown =
            ewald.computeEnergyBreakdown(config.particles, config.lattice, neighbors);

        const double realRelError = std::abs(breakdown.real - ref.eReal) / std::abs(ref.eReal);
        const double fourierRelError =
            std::abs(breakdown.reciprocal - ref.eFourier) / std::abs(ref.eFourier);
        const double selfRelError = std::abs(breakdown.self - ref.eSelf) / std::abs(ref.eSelf);
        const double intraRelError =
            std::abs(breakdown.exclusion - ref.eIntra) / std::abs(ref.eIntra);

        const double eCoulombNist = ref.eReal + ref.eFourier + ref.eSelf + ref.eIntra;
        const double eCoulombComputed =
            breakdown.real + breakdown.reciprocal + breakdown.self + breakdown.exclusion;
        const double totalRelError =
            std::abs(eCoulombComputed - eCoulombNist) / std::abs(eCoulombNist);

        // Printed unconditionally (not just on failure) -- per this
        // session's established practice, real computed-vs-reference
        // numbers must be visible regardless of pass/fail.
        std::printf(
            "%s (M=%zu, alpha=%.6f):\n"
            "  E_real:    computed=%14.4f  NIST=%14.4f  relErr=%.3e\n"
            "  E_fourier: computed=%14.4f  NIST=%14.4f  relErr=%.3e\n"
            "  E_self:    computed=%14.4f  NIST=%14.4f  relErr=%.3e\n"
            "  E_intra:   computed=%14.4f  NIST=%14.4f  relErr=%.3e\n"
            "  E_coulomb: computed=%14.4f  NIST=%14.4f  relErr=%.3e\n",
            ref.configFile.c_str(), ref.moleculeCount, alpha, breakdown.real, ref.eReal,
            realRelError, breakdown.reciprocal, ref.eFourier, fourierRelError, breakdown.self,
            ref.eSelf, selfRelError, breakdown.exclusion, ref.eIntra, intraRelError,
            eCoulombComputed, eCoulombNist, totalRelError);

        INFO("config: " << ref.configFile);
        CHECK(realRelError < kRelativeTolerance);
        CHECK(fourierRelError < kRelativeTolerance);
        CHECK(selfRelError < kRelativeTolerance);
        CHECK(intraRelError < kRelativeTolerance);
        CHECK(totalRelError < kRelativeTolerance);
    }
}
