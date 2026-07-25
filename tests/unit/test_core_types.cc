#include <catch2/catch_test_macros.hpp>

#include <array>
#include <cmath>
#include <cstdint>
#include <string>

#include "core/exceptions.hpp"
#include "core/geometry/lattice.hpp"
#include "core/math/counter_based_rng.hpp"
#include "core/math/particle_data.hpp"
#include "core/memory/aligned_allocator.hpp"
#include "core/memory/arena.hpp"

TEST_CASE("NotImplemented carries its message", "[unit][core]") {
    try {
        throw aleator::NotImplemented("Foo::bar");
    } catch (const aleator::NotImplemented& err) {
        REQUIRE(std::string(err.what()).find("Foo::bar") != std::string::npos);
    }
}

TEST_CASE("ParticleData is a plain SoA container", "[unit][core]") {
    aleator::core::ParticleData particles;
    REQUIRE(particles.empty());

    particles.push_back(1.0, 2.0, 3.0, 0.0, 0.0, 0.0, 12.0, -1.0, 0);
    REQUIRE(particles.size() == 1);
    REQUIRE(particles.x[0] == 1.0);
    REQUIRE(particles.charge[0] == -1.0);

    particles.resize(5);
    REQUIRE(particles.size() == 5);
    REQUIRE(particles.y.size() == 5);
    REQUIRE(particles.species.size() == 5);
}

TEST_CASE("Lattice stores the general triclinic matrix, not orthogonal lengths",
          "[unit][core]") {
    // A deliberately non-orthogonal (triclinic) cell.
    std::array<std::array<double, 3>, 3> matrix{};
    matrix[0] = {10.0, 0.0, 0.0};
    matrix[1] = {2.0, 8.0, 0.0};
    matrix[2] = {1.0, 1.0, 9.0};
    const aleator::core::Lattice lattice(matrix);
    // Volume = a . (b x c) for this matrix = 10 * (8*9 - 0*1) = 720.
    REQUIRE(std::abs(lattice.volume() - 720.0) < 1e-9);
}

TEST_CASE("Philox4x32Rng is declared but not yet implemented", "[unit][core]") {
    aleator::core::Philox4x32Rng rng;
    REQUIRE_THROWS_AS(rng.seed(1, 0), aleator::NotImplemented);
    REQUIRE_THROWS_AS(rng.nextU64(), aleator::NotImplemented);
    REQUIRE_THROWS_AS(rng.nextUniform(), aleator::NotImplemented);
}

TEST_CASE("AlignedAllocator returns memory aligned to the requested boundary", "[unit][core]") {
    aleator::core::AlignedAllocator<double, 64> allocator;
    double* ptr = allocator.allocate(16);
    REQUIRE(reinterpret_cast<std::uintptr_t>(ptr) % 64 == 0);
    allocator.deallocate(ptr, 16);
}

TEST_CASE("Arena hands out non-overlapping, correctly aligned allocations", "[unit][core]") {
    aleator::core::Arena arena(1024, 64);
    std::byte* a = arena.allocate(100);
    std::byte* b = arena.allocate(100);
    REQUIRE(a != nullptr);
    REQUIRE(b != nullptr);
    REQUIRE(b >= a + 100);
    REQUIRE(reinterpret_cast<std::uintptr_t>(a) % 64 == 0);

    arena.reset();
    REQUIRE(arena.used() == 0);
}
