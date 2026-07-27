// Validation anchor for core/neighbor (CLAUDE.md milestone: "assert the
// neighbor list contains EXACTLY the same pairs as brute-force O(N^2)
// within cutoff. Cubic and triclinic. Test the rebuild trigger by moving
// particles and asserting no interaction is ever missed.")
//
// Ground truth throughout is a brute-force O(N^2) minimum-image distance
// scan, implemented directly against Lattice (not reusing any part of
// CellList/VerletList), so a shared bug can't hide.

#include <catch2/catch_test_macros.hpp>

#include <algorithm>
#include <array>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <random>
#include <set>
#include <utility>
#include <vector>

#include "core/geometry/lattice.hpp"
#include "core/math/particle_data.hpp"
#include "core/neighbor/verlet_list.hpp"
#include "io/structure_io.hpp"

namespace {

using aleator::core::Lattice;
using aleator::core::ParticleData;
using aleator::core::VerletList;
using PairSet = std::set<std::pair<std::size_t, std::size_t>>;

Lattice makeLattice(std::array<double, 3> a0, std::array<double, 3> a1, std::array<double, 3> a2) {
    std::array<std::array<double, 3>, 3> m{};
    m[0] = a0;
    m[1] = a1;
    m[2] = a2;
    return Lattice(m);
}

ParticleData randomParticlesInFractionalBox(std::size_t count, unsigned seed) {
    std::mt19937 rng(seed);
    std::uniform_real_distribution<double> unit(0.0, 1.0);
    ParticleData particles;
    particles.resize(count);
    for (std::size_t i = 0; i < count; ++i) {
        particles.x[i] = unit(rng);
        particles.y[i] = unit(rng);
        particles.z[i] = unit(rng);
    }
    return particles;
}

void toCartesian(ParticleData& particles, const Lattice& lattice) {
    for (std::size_t i = 0; i < particles.size(); ++i) {
        const auto cart =
            lattice.fractionalToCartesian({particles.x[i], particles.y[i], particles.z[i]});
        particles.x[i] = cart[0];
        particles.y[i] = cart[1];
        particles.z[i] = cart[2];
    }
}

double distance(const Lattice& lattice, const ParticleData& particles, std::size_t i, std::size_t j) {
    const auto d = lattice.minimumImageDisplacement({particles.x[i], particles.y[i], particles.z[i]},
                                                      {particles.x[j], particles.y[j], particles.z[j]});
    return std::sqrt(d[0] * d[0] + d[1] * d[1] + d[2] * d[2]);
}

PairSet bruteForcePairs(const Lattice& lattice, const ParticleData& particles, double cutoff) {
    PairSet result;
    for (std::size_t i = 0; i < particles.size(); ++i) {
        for (std::size_t j = i + 1; j < particles.size(); ++j) {
            if (distance(lattice, particles, i, j) <= cutoff) {
                result.emplace(i, j);
            }
        }
    }
    return result;
}

PairSet toPairSet(const std::vector<aleator::core::NeighborPair>& pairs) {
    PairSet result;
    for (const auto& p : pairs) {
        result.emplace(p.i, p.j);
    }
    return result;
}

} // namespace

TEST_CASE("VerletList::pairs() exactly matches brute force on a random cubic system",
          "[validation][neighbor]") {
    const Lattice lattice = makeLattice({30.0, 0.0, 0.0}, {0.0, 30.0, 0.0}, {0.0, 0.0, 30.0});
    ParticleData particles = randomParticlesInFractionalBox(500, 2024);
    toCartesian(particles, lattice);

    const double cutoff = 3.0;
    VerletList list(/*skinDistance=*/1.0);
    list.build(particles, lattice, cutoff);

    REQUIRE(toPairSet(list.pairs()) == bruteForcePairs(lattice, particles, cutoff));
}

TEST_CASE("VerletList::pairs() exactly matches brute force on a random triclinic system",
          "[validation][neighbor]") {
    const Lattice lattice = makeLattice({28.0, 0.0, 0.0}, {6.0, 24.0, 0.0}, {4.0, 5.0, 22.0});
    ParticleData particles = randomParticlesInFractionalBox(500, 777);
    toCartesian(particles, lattice);

    const double cutoff = 3.0;
    VerletList list(/*skinDistance=*/1.2);
    list.build(particles, lattice, cutoff);

    REQUIRE(toPairSet(list.pairs()) == bruteForcePairs(lattice, particles, cutoff));
}

TEST_CASE("VerletList::pairs() exactly matches brute force on a real (triclinic) zeolite structure",
          "[validation][neighbor]") {
    // Reuses the real IZA PTY CIF from the io/ milestone.
    const std::string ptyCif = R"CIF(data_PTY
_cell_length_a                   7.1158
_cell_length_b                   9.1748
_cell_length_c                   9.1160
_cell_angle_alpha               84.6450
_cell_angle_beta                83.8160
_cell_angle_gamma               86.6590

loop_
_symmetry_equiv_pos_as_xyz
'+x,+y,+z'
'-x,-y,-z'

loop_
_atom_site_label
_atom_site_type_symbol
_atom_site_fract_x
_atom_site_fract_y
_atom_site_fract_z
T1    Si    0.6596    0.3547    0.1505
T2    Si    0.2383    0.3911    0.0877
T3    Si    0.0289    0.8245    0.0166
T4    Si    0.0387    0.6195    0.2988
T5    Si    0.2851    0.5520    0.5500
O1    O     0.7179    0.3434 	0.3172
O2    O     0.7193    0.5098 	0.0663
O3    O     0.7674    0.2237 	0.0642
O4    O     0.1040    0.4831 	0.2032
O5    O     0.1345    0.2467 	0.0586
O6    O     0.0185    0.7665 	0.1894
O7    O     0.1626    0.4120 	0.6082
O8    O     0.1928    0.6392 	0.4115
O9    O     0.4343    0.3421 	0.1556
O10   O     0.0000    0.0000    0.0000
O11   O     0.5000    0.5000    0.5000
)CIF";
    const auto path = std::filesystem::temp_directory_path() / "aleator_test_neighbor_PTY.cif";
    {
        std::ofstream out(path);
        out << ptyCif;
    }
    const aleator::io::StructureData structure = aleator::io::readCif(path);
    std::filesystem::remove(path);

    REQUIRE(structure.particles.size() == 30);

    // Si-O covalent bond length is ~1.6 A; O...O across a tetrahedron is
    // ~2.6 A. A 3.2 A cutoff is well within L_perp/2 for this cell
    // (shortest axis ~7.1 A) and gives a non-trivial pair count.
    const double cutoff = 3.2;
    VerletList list(/*skinDistance=*/1.0);
    list.build(structure.particles, structure.lattice, cutoff);

    const auto brute = bruteForcePairs(structure.lattice, structure.particles, cutoff);
    REQUIRE_FALSE(brute.empty());
    REQUIRE(toPairSet(list.pairs()) == brute);
}

TEST_CASE(
    "VerletList's displacement-based rebuild trigger keeps pairs() exact through continuous "
    "random motion without missing any interaction",
    "[validation][neighbor]") {
    const Lattice lattice = makeLattice({25.0, 0.0, 0.0}, {5.0, 22.0, 0.0}, {3.0, 4.0, 20.0});
    // Kept deliberately small (not for algorithmic reasons — the brute-force
    // ground truth this test recomputes every single iteration is O(N^2),
    // and this test calls it ~100x, on top of exercising
    // minimumImageDisplacement millions of times total in an unoptimized
    // Debug build. 60 particles x 120 iterations is still a real,
    // exact-every-step check of the rebuild trigger; it just doesn't need
    // to be 200 x 500 to prove that.
    ParticleData particles = randomParticlesInFractionalBox(60, 555);
    toCartesian(particles, lattice);

    const double cutoff = 3.0;
    const double skin = 1.0;
    VerletList list(skin);
    list.build(particles, lattice, cutoff);

    std::mt19937 rng(9001);
    // Steps deliberately small relative to the skin, so most steps do NOT
    // require a rebuild — exercising the "stay valid between rebuilds"
    // path, not just "rebuild every step".
    std::uniform_real_distribution<double> step(-0.05, 0.05);

    int rebuildCount = 0;
    for (int iteration = 0; iteration < 120; ++iteration) {
        for (std::size_t i = 0; i < particles.size(); ++i) {
            particles.x[i] += step(rng);
            particles.y[i] += step(rng);
            particles.z[i] += step(rng);
        }

        if (list.needsRebuild(particles)) {
            list.build(particles, lattice, cutoff);
            ++rebuildCount;
        }

        const auto brute = bruteForcePairs(lattice, particles, cutoff);
        const auto tracked = toPairSet(list.pairs());
        INFO("iteration " << iteration << ", rebuilds so far " << rebuildCount);
        REQUIRE(tracked == brute);
    }

    // The trigger should have fired at least once (otherwise this test
    // isn't exercising the mechanism at all) but not every single
    // iteration (otherwise the skin buffer isn't doing anything either).
    REQUIRE(rebuildCount > 0);
    REQUIRE(rebuildCount < 120);
}

TEST_CASE("Ignoring needsRebuild() can genuinely miss an interaction — the trigger is necessary, "
          "not decorative",
          "[validation][neighbor]") {
    // Two particles placed just outside cutoff + skin at build time, in an
    // otherwise-empty system, then one is moved a large distance straight
    // at the other without ever rebuilding. needsRebuild() must catch
    // this; and if the caller perversely ignores it, pairs() must indeed
    // be missing the interaction that has since formed, proving the skin
    // mechanism is load-bearing rather than a no-op.
    const Lattice lattice = makeLattice({50.0, 0.0, 0.0}, {0.0, 50.0, 0.0}, {0.0, 0.0, 50.0});
    const double cutoff = 3.0;
    const double skin = 1.0;

    ParticleData particles;
    particles.resize(2);
    particles.x[0] = 10.0;
    particles.y[0] = 10.0;
    particles.z[0] = 10.0;
    particles.x[1] = 15.0; // 5.0 A away: outside cutoff+skin (4.0)
    particles.y[1] = 10.0;
    particles.z[1] = 10.0;

    VerletList list(skin);
    list.build(particles, lattice, cutoff);
    REQUIRE(list.pairs().empty());

    // Move particle 1 well past what a half-skin-per-particle budget
    // allows in one shot, landing it within cutoff of particle 0.
    particles.x[1] = 11.5; // now 1.5 A away: within cutoff

    REQUIRE(list.needsRebuild(particles)); // the trigger correctly fires

    // If the caller ignores that and asks for pairs() anyway:
    REQUIRE(list.pairs().empty()); // ...the interaction really is missing.

    // Honoring the trigger fixes it.
    list.build(particles, lattice, cutoff);
    REQUIRE(list.pairs().size() == 1);
}
