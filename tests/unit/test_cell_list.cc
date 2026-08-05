// core/neighbor milestone: CellList is a candidate generator, not a final
// answer — these tests check its two required properties directly:
// (1) coverage: every truly-within-radius pair (by an independent
//     brute-force minimum-image check) appears among the candidates;
// (2) no duplicates: each unordered pair is visited at most once,
//     including at small grid dimensions (nx/ny/nz == 1 or 2), which is
//     exactly where a careless wraparound offset scheme double-visits a
//     cell pair (see the comment in cell_list.cc).

#include <catch2/catch_test_macros.hpp>

#include <algorithm>
#include <array>
#include <cmath>
#include <random>
#include <set>
#include <stdexcept>
#include <utility>
#include <vector>

#include "core/neighbor/cell_list.hpp"

namespace {

using aleator::core::CellList;
using aleator::core::Lattice;
using aleator::core::ParticleData;
using PairSet = std::set<std::pair<std::size_t, std::size_t>>;

Lattice makeLattice(std::array<double, 3> a0, std::array<double, 3> a1, std::array<double, 3> a2) {
    std::array<std::array<double, 3>, 3> m{};
    m[0] = a0;
    m[1] = a1;
    m[2] = a2;
    return Lattice(m);
}

ParticleData randomParticles(std::size_t count, unsigned seed) {
    std::mt19937 rng(seed);
    std::uniform_real_distribution<double> unit(0.0, 1.0);
    ParticleData particles;
    particles.resize(count);
    for (std::size_t i = 0; i < count; ++i) {
        particles.x[i] = unit(rng);
        particles.y[i] = unit(rng);
        particles.z[i] = unit(rng);
        particles.species[i] = 0;
    }
    return particles;
}

/// Converts fractional-space random particles (in [0,1)^3, as produced by
/// randomParticles) to Cartesian under `lattice`, in place.
void toCartesian(ParticleData& particles, const Lattice& lattice) {
    for (std::size_t i = 0; i < particles.size(); ++i) {
        const auto cart = lattice.fractionalToCartesian({particles.x[i], particles.y[i], particles.z[i]});
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

PairSet bruteForceWithinRadius(const Lattice& lattice, const ParticleData& particles, double radius) {
    PairSet result;
    for (std::size_t i = 0; i < particles.size(); ++i) {
        for (std::size_t j = i + 1; j < particles.size(); ++j) {
            if (distance(lattice, particles, i, j) <= radius) {
                result.emplace(i, j);
            }
        }
    }
    return result;
}

/// Runs the coverage + no-duplicates checks for one (lattice, radius,
/// particle count) configuration and returns the grid dimensions actually
/// used, so callers can assert they hit the specific grid size they meant
/// to exercise.
std::array<int, 3> checkCellList(const Lattice& lattice, const ParticleData& particles,
                                  double radius) {
    CellList cellList;
    cellList.build(particles, lattice, radius);

    std::vector<std::pair<std::size_t, std::size_t>> candidates;
    cellList.forEachCandidatePair([&](std::size_t i, std::size_t j) {
        REQUIRE(i < j);
        candidates.emplace_back(i, j);
    });

    const PairSet candidateSet(candidates.begin(), candidates.end());
    REQUIRE(candidateSet.size() == candidates.size()); // no duplicate candidates

    const PairSet truePairs = bruteForceWithinRadius(lattice, particles, radius);
    for (const auto& pair : truePairs) {
        INFO("missing candidate pair (" << pair.first << ", " << pair.second << ")");
        REQUIRE(candidateSet.count(pair) == 1);
    }

    return {cellList.cellCountX(), cellList.cellCountY(), cellList.cellCountZ()};
}

} // namespace

TEST_CASE("CellList finds every within-radius pair with no duplicates on a cubic cell",
          "[unit][neighbor]") {
    const Lattice lattice = makeLattice({20.0, 0.0, 0.0}, {0.0, 20.0, 0.0}, {0.0, 0.0, 20.0});
    ParticleData particles = randomParticles(300, 42);
    toCartesian(particles, lattice);

    // radius = 2.0 against a 20 A box gives a comfortably large (~10^3) grid.
    const auto grid = checkCellList(lattice, particles, 2.0);
    REQUIRE(grid[0] >= 3);
}

TEST_CASE("CellList is correct when forced to a single cell along every axis (nx=ny=nz=1)",
          "[unit][neighbor]") {
    const Lattice lattice = makeLattice({10.0, 0.0, 0.0}, {0.0, 10.0, 0.0}, {0.0, 0.0, 10.0});
    ParticleData particles = randomParticles(80, 7);
    toCartesian(particles, lattice);

    // radius larger than the whole cell's perpendicular width forces n=1.
    const auto grid = checkCellList(lattice, particles, 12.0);
    REQUIRE(grid[0] == 1);
    REQUIRE(grid[1] == 1);
    REQUIRE(grid[2] == 1);
}

TEST_CASE("CellList is correct when forced to exactly 2 cells along every axis",
          "[unit][neighbor]") {
    const Lattice lattice = makeLattice({10.0, 0.0, 0.0}, {0.0, 10.0, 0.0}, {0.0, 0.0, 10.0});
    ParticleData particles = randomParticles(150, 99);
    toCartesian(particles, lattice);

    // perpendicularWidth / radius in [2,3) forces n=2 along every axis:
    // 10 / 4.0 = 2.5 -> floor = 2.
    const auto grid = checkCellList(lattice, particles, 4.0);
    REQUIRE(grid[0] == 2);
    REQUIRE(grid[1] == 2);
    REQUIRE(grid[2] == 2);
}

TEST_CASE("CellList is correct on a triclinic cell", "[unit][neighbor]") {
    const Lattice lattice = makeLattice({15.0, 0.0, 0.0}, {4.0, 12.0, 0.0}, {2.5, 3.0, 11.0});
    ParticleData particles = randomParticles(300, 123);
    toCartesian(particles, lattice);

    checkCellList(lattice, particles, 2.5);
}

TEST_CASE("CellList handles zero and one particle without crashing", "[unit][neighbor]") {
    const Lattice lattice = makeLattice({10.0, 0.0, 0.0}, {0.0, 10.0, 0.0}, {0.0, 0.0, 10.0});

    ParticleData empty;
    CellList emptyList;
    REQUIRE_NOTHROW(emptyList.build(empty, lattice, 2.0));
    int visits = 0;
    emptyList.forEachCandidatePair([&](std::size_t, std::size_t) { ++visits; });
    REQUIRE(visits == 0);

    ParticleData single;
    single.resize(1);
    single.x[0] = 1.0;
    single.y[0] = 1.0;
    single.z[0] = 1.0;
    CellList singleList;
    REQUIRE_NOTHROW(singleList.build(single, lattice, 2.0));
    visits = 0;
    singleList.forEachCandidatePair([&](std::size_t, std::size_t) { ++visits; });
    REQUIRE(visits == 0);
}

namespace {

/// Brute-force set of indices within `radius` of an arbitrary (not
/// necessarily binned) query point, the ground truth for
/// forEachIndexNear()'s coverage guarantee.
std::set<std::size_t> bruteForceIndicesNear(const Lattice& lattice, const ParticleData& particles,
                                             const std::array<double, 3>& queryPoint,
                                             double radius) {
    std::set<std::size_t> result;
    for (std::size_t i = 0; i < particles.size(); ++i) {
        const auto d = lattice.minimumImageDisplacement(
            queryPoint, {particles.x[i], particles.y[i], particles.z[i]});
        const double r = std::sqrt(d[0] * d[0] + d[1] * d[1] + d[2] * d[2]);
        if (r <= radius) {
            result.insert(i);
        }
    }
    return result;
}

} // namespace

TEST_CASE("CellList::forEachIndexNear finds every within-radius index for arbitrary query points",
          "[unit][neighbor]") {
    const Lattice lattice = makeLattice({20.0, 0.0, 0.0}, {0.0, 20.0, 0.0}, {0.0, 0.0, 20.0});
    ParticleData particles = randomParticles(300, 55);
    toCartesian(particles, lattice);

    const double radius = 2.0;
    CellList cellList;
    cellList.build(particles, lattice, radius);

    std::mt19937 queryRng(4242);
    std::uniform_real_distribution<double> unit(0.0, 1.0);
    for (int trial = 0; trial < 50; ++trial) {
        const auto queryPoint = lattice.fractionalToCartesian(
            {unit(queryRng), unit(queryRng), unit(queryRng)});

        std::vector<std::size_t> candidates;
        cellList.forEachIndexNear(lattice, queryPoint,
                                   [&](std::size_t j) { candidates.push_back(j); });
        const std::set<std::size_t> candidateSet(candidates.begin(), candidates.end());
        REQUIRE(candidateSet.size() == candidates.size()); // no duplicates

        const auto truth = bruteForceIndicesNear(lattice, particles, queryPoint, radius);
        for (std::size_t i : truth) {
            INFO("missing candidate index " << i << " for trial " << trial);
            REQUIRE(candidateSet.count(i) == 1);
        }
    }
}

TEST_CASE("CellList::forEachIndexNear is correct on a triclinic cell", "[unit][neighbor]") {
    const Lattice lattice = makeLattice({15.0, 0.0, 0.0}, {4.0, 12.0, 0.0}, {2.5, 3.0, 11.0});
    ParticleData particles = randomParticles(300, 321);
    toCartesian(particles, lattice);

    const double radius = 2.5;
    CellList cellList;
    cellList.build(particles, lattice, radius);

    std::mt19937 queryRng(99);
    std::uniform_real_distribution<double> unit(0.0, 1.0);
    for (int trial = 0; trial < 50; ++trial) {
        const auto queryPoint =
            lattice.fractionalToCartesian({unit(queryRng), unit(queryRng), unit(queryRng)});
        std::vector<std::size_t> candidates;
        cellList.forEachIndexNear(lattice, queryPoint,
                                   [&](std::size_t j) { candidates.push_back(j); });
        const std::set<std::size_t> candidateSet(candidates.begin(), candidates.end());

        const auto truth = bruteForceIndicesNear(lattice, particles, queryPoint, radius);
        for (std::size_t i : truth) {
            INFO("missing candidate index " << i << " for trial " << trial);
            REQUIRE(candidateSet.count(i) == 1);
        }
    }
}

TEST_CASE("CellList::build with startIndex only bins particles at or above startIndex",
          "[unit][neighbor]") {
    const Lattice lattice = makeLattice({20.0, 0.0, 0.0}, {0.0, 20.0, 0.0}, {0.0, 0.0, 20.0});
    ParticleData particles = randomParticles(300, 7);
    toCartesian(particles, lattice);

    constexpr std::size_t kStartIndex = 200;
    CellList cellList;
    cellList.build(particles, lattice, 2.0, kStartIndex);

    for (std::size_t i = 0; i < particles.size(); ++i) {
        std::vector<std::size_t> candidates;
        cellList.forEachIndexNear(lattice, {particles.x[i], particles.y[i], particles.z[i]},
                                   [&](std::size_t j) { candidates.push_back(j); });
        for (std::size_t j : candidates) {
            REQUIRE(j >= kStartIndex);
        }
    }
    // The particle's own bin, queried at its own position, must include
    // itself whenever its index is >= kStartIndex (sanity check that
    // startIndex didn't drop everything).
    std::vector<std::size_t> selfCandidates;
    cellList.forEachIndexNear(
        lattice, {particles.x[250], particles.y[250], particles.z[250]},
        [&](std::size_t j) { selfCandidates.push_back(j); });
    REQUIRE(std::find(selfCandidates.begin(), selfCandidates.end(), std::size_t{250}) !=
            selfCandidates.end());
}

TEST_CASE("CellList::build rejects a non-positive radius", "[unit][neighbor]") {
    const Lattice lattice = makeLattice({10.0, 0.0, 0.0}, {0.0, 10.0, 0.0}, {0.0, 0.0, 10.0});
    ParticleData particles = randomParticles(10, 1);
    CellList cellList;
    REQUIRE_THROWS_AS(cellList.build(particles, lattice, 0.0), std::invalid_argument);
    REQUIRE_THROWS_AS(cellList.build(particles, lattice, -1.0), std::invalid_argument);
}
