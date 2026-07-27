// Validation anchor (CLAUDE.md #4: "LJ energy + cutoff | NIST SRSW LJ
// reference configurations | 1e-8 relative"). Uses the real published NIST
// SRSW Lennard-Jones fluid reference data — sample configurations fetched
// from nist.gov (see data/nist_lj/PROVENANCE.md; not invented) and the
// reference U_pair / W_pair / U_LRC values transcribed from the
// "Reference Calculations for the Internal Energy and Pair Virial" table
// at https://www.nist.gov/mml/csd/chemical-informatics-group/lennard-jones-fluid-reference-calculations-cuboid-cell.
//
// A note on the "1e-8 relative" figure in CLAUDE.md's validation-anchor
// table: NIST's published reference values are given to 5 significant
// figures (e.g. "-4.3515E+03"), so 1e-8 relative agreement against them is
// not something any implementation could demonstrate — the reference
// itself only carries ~1e-5 relative precision (rounding in the last
// digit). What IS meaningful, and what this test actually checks, is
// agreement to within that ~1e-5 rounding budget (tolerance set to 5e-5,
// slightly looser to avoid spurious failures from worst-case rounding);
// double-precision arithmetic throughout means the code's own internal
// precision is far tighter than that, but there is nothing to compare it
// against beyond the published digits. Real measured relative errors are
// printed in each failure message via INFO/CAPTURE so this is auditable.

#include <catch2/catch_test_macros.hpp>

#include <array>
#include <cmath>
#include <cstdint>
#include <filesystem>
#include <fstream>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include "core/geometry/lattice.hpp"
#include "core/math/particle_data.hpp"
#include "core/neighbor/verlet_list.hpp"
#include "forcefield/pairwise/lennard_jones.hpp"

namespace {

using aleator::core::Lattice;
using aleator::core::ParticleData;
using aleator::core::VerletList;
using aleator::forcefield::LennardJones;
using aleator::forcefield::LennardJonesParameters;
using aleator::forcefield::LennardJonesTruncation;

// Loosened from the reference data's own ~1e-5 rounding precision — see
// file header comment.
constexpr double kRelativeTolerance = 5e-5;

struct NistConfig {
    Lattice lattice;
    ParticleData particles;
};

/// Parses a NIST SRSW sample configuration file: line 1 is the box x/y/z
/// extent, line 2 is the particle count N, followed by N "id x y z" rows.
/// All values are whitespace-separated; scientific notation ("1.0E+01")
/// parses directly via operator>>.
NistConfig loadNistConfig(const std::filesystem::path& path) {
    std::ifstream in(path);
    if (!in) {
        throw std::runtime_error("could not open NIST config file: " + path.string());
    }

    double lx = 0.0;
    double ly = 0.0;
    double lz = 0.0;
    in >> lx >> ly >> lz;

    std::size_t n = 0;
    in >> n;

    ParticleData particles;
    particles.resize(n);
    for (std::size_t k = 0; k < n; ++k) {
        std::size_t id = 0;
        double x = 0.0;
        double y = 0.0;
        double z = 0.0;
        in >> id >> x >> y >> z;
        if (!in) {
            throw std::runtime_error("malformed NIST config file (short read) at particle " +
                                      std::to_string(k) + ": " + path.string());
        }
        particles.x[k] = x;
        particles.y[k] = y;
        particles.z[k] = z;
        particles.species[k] = 0;
    }

    std::array<std::array<double, 3>, 3> m{};
    m[0] = {lx, 0.0, 0.0};
    m[1] = {0.0, ly, 0.0};
    m[2] = {0.0, 0.0, lz};
    return {Lattice(m), std::move(particles)};
}

std::filesystem::path nistDataFile(const std::string& name) {
    return std::filesystem::path(ALEATOR_TEST_DATA_DIR) / "data" / "nist_lj" / name;
}

struct SchemeReference {
    std::string label;
    double cutoff;
    LennardJonesTruncation truncation;
    double upair;
    double wpair;
    double ulrc; // only meaningful (and only checked) when addTailCorrection is true
    bool addTailCorrection;
};

struct ConfigReference {
    std::string configFile;
    std::vector<SchemeReference> schemes;
};

// Transcribed from the NIST reference table (see file header for the exact
// source URL). epsilon = sigma = 1 (LJ reduced units) for every particle in
// every configuration.
const std::vector<ConfigReference>& nistReferenceData() {
    static const std::vector<ConfigReference> data{
        {"lj_sample_config_periodic1.txt",
         {
             {"LRC rc*=3.0", 3.0, LennardJonesTruncation::Truncated, -4.3515e+03, -5.6867e+02,
              -1.9849e+02, true},
             {"LRC rc*=4.0", 4.0, LennardJonesTruncation::Truncated, -4.4675e+03, -1.2639e+03,
              -8.3769e+01, true},
             {"LFS rc*=3.0", 3.0, LennardJonesTruncation::LinearForceShifted, -3.8709e+03,
              3.1754e+02, 0.0, false},
         }},
        {"lj_sample_config_periodic2.txt",
         {
             {"LRC rc*=3.0", 3.0, LennardJonesTruncation::Truncated, -6.9000e+02, -5.6846e+02,
              -2.4230e+01, true},
             {"LRC rc*=4.0", 4.0, LennardJonesTruncation::Truncated, -7.0460e+02, -6.5599e+02,
              -1.0226e+01, true},
             {"LFS rc*=3.0", 3.0, LennardJonesTruncation::LinearForceShifted, -6.2012e+02,
              -4.4533e+02, 0.0, false},
         }},
        {"lj_sample_config_periodic3.txt",
         {
             {"LRC rc*=3.0", 3.0, LennardJonesTruncation::Truncated, -1.1467e+03, -1.1649e+03,
              -4.9622e+01, true},
             {"LRC rc*=4.0", 4.0, LennardJonesTruncation::Truncated, -1.1754e+03, -1.3371e+03,
              -2.0942e+01, true},
             {"LFS rc*=3.0", 3.0, LennardJonesTruncation::LinearForceShifted, -1.0210e+03,
              -9.3578e+02, 0.0, false},
         }},
        {"lj_sample_config_periodic4.txt",
         {
             {"LRC rc*=3.0", 3.0, LennardJonesTruncation::Truncated, -1.6790e+01, -4.6249e+01,
              -5.4517e-01, true},
             {"LRC rc*=4.0", 4.0, LennardJonesTruncation::Truncated, -1.7060e+01, -4.7869e+01,
              -2.3008e-01, true},
             {"LFS rc*=3.0", 3.0, LennardJonesTruncation::LinearForceShifted, -1.5001e+01,
              -4.3096e+01, 0.0, false},
         }},
    };
    return data;
}

} // namespace

TEST_CASE("LennardJones matches the published NIST SRSW reference energies and virials",
          "[validation][forcefield]") {
    for (const auto& configRef : nistReferenceData()) {
        const NistConfig config = loadNistConfig(nistDataFile(configRef.configFile));

        for (const auto& scheme : configRef.schemes) {
            INFO("config file: " << configRef.configFile << ", scheme: " << scheme.label
                                  << ", N=" << config.particles.size());

            const std::vector<LennardJonesParameters> species{{1.0, 1.0}};
            const LennardJones lj(species, scheme.cutoff, scheme.truncation);

            VerletList neighbors(/*skinDistance=*/0.5);
            neighbors.build(config.particles, config.lattice, scheme.cutoff);

            const double upair = lj.computeEnergy(config.particles, config.lattice, neighbors);
            const double wpair = lj.computeVirial(config.particles, config.lattice, neighbors);

            const double upairRelError = std::abs(upair - scheme.upair) / std::abs(scheme.upair);
            const double wpairRelError = std::abs(wpair - scheme.wpair) / std::abs(scheme.wpair);
            INFO("U_pair: computed=" << upair << " NIST=" << scheme.upair
                                      << " relError=" << upairRelError);
            INFO("W_pair: computed=" << wpair << " NIST=" << scheme.wpair
                                      << " relError=" << wpairRelError);
            REQUIRE(upairRelError < kRelativeTolerance);
            REQUIRE(wpairRelError < kRelativeTolerance);

            if (scheme.addTailCorrection) {
                const double ulrc = lj.tailEnergyCorrection(config.particles, config.lattice);
                const double ulrcRelError = std::abs(ulrc - scheme.ulrc) / std::abs(scheme.ulrc);
                INFO("U_LRC: computed=" << ulrc << " NIST=" << scheme.ulrc
                                         << " relError=" << ulrcRelError);
                REQUIRE(ulrcRelError < kRelativeTolerance);
            }
        }
    }
}
