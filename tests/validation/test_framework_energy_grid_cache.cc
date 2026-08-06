// Validation anchor for CLAUDE.md section 5's supercell-era
// FrameworkEnergyGrid milestone: "Serialize FrameworkEnergyGrid to disk...
// A stale grid silently producing wrong energies is the worst failure mode
// in this project -- a mismatched hash must be a hard error, never a
// silent rebuild under a different key, and never a silent reuse."
//
// Two things are tested, deliberately kept separate:
//   1. Round-trip: a deserialized grid gives BIT-IDENTICAL energies to the
//      in-memory one it came from -- raw bits (operator==), not a
//      tolerance, since serialize()/loadFromCache() round-trip raw stored
//      doubles with no arithmetic in between.
//   2. Poisoning: mutating any single field of the cache key changes the
//      content-addressed cache filename, so a lookup under the mutated key
//      is a normal cache MISS against the original file (never a match) --
//      and, separately, if a file that does not actually match a
//      requested key is ever encountered at that key's expected path (the
//      hash-collision/corruption/tampering scenario), loading it is a HARD
//      ERROR, never a silent reuse and never a silent "just rebuild as if
//      nothing were wrong".

#include <catch2/catch_test_macros.hpp>

#include <cstdio>
#include <filesystem>
#include <fstream>
#include <stdexcept>
#include <string>
#include <vector>

#include "core/geometry/lattice.hpp"
#include "core/math/particle_data.hpp"
#include "engines/monte_carlo/framework_energy_grid.hpp"
#include "forcefield/pairwise/lennard_jones.hpp"

using aleator::core::Lattice;
using aleator::core::ParticleData;
using aleator::engines::FrameworkEnergyGrid;
using aleator::engines::FrameworkEnergyGridCacheKey;
using aleator::engines::FrameworkEnergyGridCacheRequest;
using aleator::engines::loadOrBuildFrameworkEnergyGrid;
using aleator::engines::makeFrameworkEnergyGridCacheKey;
using aleator::forcefield::LennardJones;
using aleator::forcefield::LennardJonesParameters;
using aleator::forcefield::LennardJonesTruncation;

namespace {

/// A small, fast synthetic system -- this file is about cache-key
/// mechanics, not accuracy or a specific real framework, so keeping it
/// small keeps every test case here fast (the real-IRMOF-1 accuracy sweep
/// lives in test_framework_energy_grid.cc).
struct TestSystem {
    ParticleData particles;
    Lattice lattice;
    std::vector<LennardJonesParameters> ljParameters;
    std::uint32_t guestSpeciesIndex;
};

TestSystem buildTestSystem() {
    TestSystem sys;
    std::array<std::array<double, 3>, 3> m{};
    m[0] = {15.0, 0.0, 0.0};
    m[1] = {0.0, 15.0, 0.0};
    m[2] = {0.0, 0.0, 15.0};
    sys.lattice = Lattice(m);

    sys.particles.push_back(2.0, 2.0, 2.0, 0, 0, 0, 65.38, 0.0, 0);
    sys.particles.push_back(7.0, 7.0, 7.0, 0, 0, 0, 16.0, 0.0, 1);
    sys.particles.push_back(12.0, 3.0, 9.0, 0, 0, 0, 65.38, 0.0, 0);

    sys.ljParameters = {{62.4, 2.46}, {48.2, 3.03}}; // species 0 (Zn-like), 1 (O-like)
    sys.guestSpeciesIndex = static_cast<std::uint32_t>(sys.ljParameters.size());
    sys.ljParameters.push_back({158.5, 3.72}); // guest species
    return sys;
}

} // namespace

TEST_CASE("FrameworkEnergyGrid cache round-trip is bit-identical, not just close",
          "[validation][montecarlo][cache]") {
    const TestSystem sys = buildTestSystem();
    const LennardJones forceField(sys.ljParameters, 6.0, LennardJonesTruncation::Shifted);

    const auto cacheDir =
        std::filesystem::temp_directory_path() / "aleator_test_grid_cache_roundtrip";
    std::filesystem::remove_all(cacheDir);

    FrameworkEnergyGridCacheRequest request{cacheDir, {sys.guestSpeciesIndex}, 1.0};

    bool firstHit = true;
    const auto original = loadOrBuildFrameworkEnergyGrid(forceField, sys.particles, sys.lattice,
                                                           request, &firstHit);
    REQUIRE_FALSE(firstHit); // nothing cached yet

    bool secondHit = false;
    const auto reloaded = loadOrBuildFrameworkEnergyGrid(forceField, sys.particles, sys.lattice,
                                                           request, &secondHit);
    REQUIRE(secondHit); // must come from the file just written, not a rebuild

    // Sample many points, including off-node positions (exercising the
    // trilinear interpolation path, not just stored node values).
    int mismatches = 0;
    for (int i = 0; i < 50; ++i) {
        const std::array<double, 3> p{1.0 + i * 0.27, 2.0 + i * 0.19, 3.0 + i * 0.11};
        const double e1 = original.interpolate(sys.guestSpeciesIndex, p);
        const double e2 = reloaded.interpolate(sys.guestSpeciesIndex, p);
        if (e1 != e2) { // bit-exact, deliberately not a tolerance comparison
            ++mismatches;
            INFO("point " << i << ": original=" << e1 << " reloaded=" << e2);
        }
    }
    CHECK(mismatches == 0);

    std::filesystem::remove_all(cacheDir);
}

TEST_CASE("FrameworkEnergyGrid cache: mutating any single key field is a real cache miss",
          "[validation][montecarlo][cache]") {
    const TestSystem sys = buildTestSystem();
    const LennardJones forceField(sys.ljParameters, 6.0, LennardJonesTruncation::Shifted);

    const auto cacheDir = std::filesystem::temp_directory_path() / "aleator_test_grid_cache_poison";
    std::filesystem::remove_all(cacheDir);

    const auto baselineKey = makeFrameworkEnergyGridCacheKey(forceField, sys.particles, sys.lattice,
                                                               {sys.guestSpeciesIndex}, 1.0);
    // Populate the cache under the baseline key.
    FrameworkEnergyGridCacheRequest request{cacheDir, {sys.guestSpeciesIndex}, 1.0};
    (void)loadOrBuildFrameworkEnergyGrid(forceField, sys.particles, sys.lattice, request);
    REQUIRE(FrameworkEnergyGrid::cacheHit(cacheDir, baselineKey));

    auto expectMiss = [&](const FrameworkEnergyGridCacheKey& mutated, const std::string& label) {
        INFO("mutated field: " << label);
        CHECK_FALSE(FrameworkEnergyGrid::cacheHit(cacheDir, mutated));
        CHECK(baselineKey.hash() != mutated.hash());
    };

    {
        auto k = baselineKey;
        k.atomX[0] += 0.5;
        expectMiss(k, "atomX");
    }
    {
        auto k = baselineKey;
        k.atomY[1] -= 0.3;
        expectMiss(k, "atomY");
    }
    {
        auto k = baselineKey;
        k.atomZ[2] += 1.1;
        expectMiss(k, "atomZ");
    }
    {
        auto k = baselineKey;
        k.atomSpecies[0] = k.atomSpecies[0] == 0 ? 1 : 0;
        expectMiss(k, "atomSpecies");
    }
    {
        auto k = baselineKey;
        k.latticeMatrix[0][0] += 0.2;
        expectMiss(k, "latticeMatrix");
    }
    {
        auto k = baselineKey;
        k.perSpeciesLJ[0].epsilon *= 1.01;
        expectMiss(k, "perSpeciesLJ[].epsilon");
    }
    {
        auto k = baselineKey;
        k.perSpeciesLJ[0].sigma *= 1.01;
        expectMiss(k, "perSpeciesLJ[].sigma");
    }
    {
        auto k = baselineKey;
        k.cutoffAngstrom += 0.5;
        expectMiss(k, "cutoffAngstrom");
    }
    {
        auto k = baselineKey;
        k.truncation = static_cast<std::uint8_t>(LennardJonesTruncation::Truncated);
        REQUIRE(k.truncation != baselineKey.truncation); // sanity: baseline was Shifted
        expectMiss(k, "truncation");
    }
    {
        auto k = baselineKey;
        k.guestSpeciesIds.push_back(99); // a second, different guest species
        expectMiss(k, "guestSpeciesIds");
    }
    {
        auto k = baselineKey;
        k.spacingAngstrom += 0.25;
        expectMiss(k, "spacingAngstrom");
    }
    {
        auto k = baselineKey;
        k.energyCapKelvin *= 2.0;
        expectMiss(k, "energyCapKelvin");
    }

    std::filesystem::remove_all(cacheDir);
}

TEST_CASE("FrameworkEnergyGrid cache: a file whose stored key doesn't match the requested key is "
          "a hard error, never a silent reuse",
          "[validation][montecarlo][cache]") {
    const TestSystem sys = buildTestSystem();
    const LennardJones forceField(sys.ljParameters, 6.0, LennardJonesTruncation::Shifted);

    const auto cacheDir =
        std::filesystem::temp_directory_path() / "aleator_test_grid_cache_collision";
    std::filesystem::remove_all(cacheDir);

    const auto baselineKey = makeFrameworkEnergyGridCacheKey(forceField, sys.particles, sys.lattice,
                                                               {sys.guestSpeciesIndex}, 1.0);
    FrameworkEnergyGridCacheRequest request{cacheDir, {sys.guestSpeciesIndex}, 1.0};
    (void)loadOrBuildFrameworkEnergyGrid(forceField, sys.particles, sys.lattice, request);

    // Simulate exactly the dangerous scenario the milestone names: a file
    // exists at the path a DIFFERENT key would look up (as if a hash
    // collision, a corrupted directory, or a foreign cache had put it
    // there), and its actual stored contents don't match what's being
    // asked for. Forced here by copying the real (baseline-key) file to
    // sit at a different key's computed path, rather than relying on
    // actually finding a real FNV-1a collision.
    auto differentKey = baselineKey;
    differentKey.cutoffAngstrom += 1.0;
    const auto baselinePath = FrameworkEnergyGrid::cacheFilePath(cacheDir, baselineKey);
    const auto collidingPath = FrameworkEnergyGrid::cacheFilePath(cacheDir, differentKey);
    REQUIRE(baselinePath != collidingPath);
    std::filesystem::copy_file(baselinePath, collidingPath);

    // Looking up `differentKey` finds a file (the one just copied in), but
    // that file's stored key is `baselineKey`, not `differentKey` -- must
    // hard-error, not silently report a miss and not silently return the
    // (wrong) grid.
    CHECK_THROWS_AS(FrameworkEnergyGrid::cacheHit(cacheDir, differentKey), std::runtime_error);
    CHECK_THROWS_AS(FrameworkEnergyGrid::loadFromCache(cacheDir, differentKey), std::runtime_error);

    std::filesystem::remove_all(cacheDir);
}

TEST_CASE("FrameworkEnergyGrid cache: an on-disk format version mismatch is a hard error",
          "[validation][montecarlo][cache]") {
    const TestSystem sys = buildTestSystem();
    const LennardJones forceField(sys.ljParameters, 6.0, LennardJonesTruncation::Shifted);

    const auto cacheDir =
        std::filesystem::temp_directory_path() / "aleator_test_grid_cache_version";
    std::filesystem::remove_all(cacheDir);

    const auto key = makeFrameworkEnergyGridCacheKey(forceField, sys.particles, sys.lattice,
                                                       {sys.guestSpeciesIndex}, 1.0);
    FrameworkEnergyGridCacheRequest request{cacheDir, {sys.guestSpeciesIndex}, 1.0};
    (void)loadOrBuildFrameworkEnergyGrid(forceField, sys.particles, sys.lattice, request);

    const auto path = FrameworkEnergyGrid::cacheFilePath(cacheDir, key);
    {
        // Overwrite just the leading 4-byte format-version field in place
        // with a value that can never be the real (current) version.
        std::fstream file(path, std::ios::binary | std::ios::in | std::ios::out);
        REQUIRE(file.is_open());
        const std::uint32_t bogusVersion = 0xffffffffu;
        file.write(reinterpret_cast<const char*>(&bogusVersion), sizeof(bogusVersion));
    }

    CHECK_THROWS_AS(FrameworkEnergyGrid::cacheHit(cacheDir, key), std::runtime_error);
    CHECK_THROWS_AS(FrameworkEnergyGrid::loadFromCache(cacheDir, key), std::runtime_error);

    std::filesystem::remove_all(cacheDir);
}
