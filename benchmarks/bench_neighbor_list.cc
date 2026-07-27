// core/neighbor milestone: honest scalar baselines for CellList/VerletList
// build and traversal at N = 1e3..1e6. Deliberately NOT optimized —
// CLAUDE.md #5: profile first, optimize only where the profiler points,
// and this is the "record the baseline" step that has to happen before
// that. Run via the `bench` preset for a realistic (Release,
// -march=native, scoped to this binary only) number; a Debug build number
// here would not be representative of anything.
//
// Density is held constant across N (box volume scales with N) so the
// average neighbor count per particle — and thus the per-particle work —
// stays comparable across the whole N range; otherwise growing N at fixed
// box size would silently also be growing the density, confounding the
// scaling picture.

#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <random>

#include <benchmark/benchmark.h>

#include "core/geometry/lattice.hpp"
#include "core/math/particle_data.hpp"
#include "core/neighbor/cell_list.hpp"
#include "core/neighbor/verlet_list.hpp"

namespace {

using aleator::core::CellList;
using aleator::core::Lattice;
using aleator::core::ParticleData;
using aleator::core::VerletList;

constexpr double kCutoff = 10.0;      // A, a representative LJ-scale cutoff
constexpr double kSkin = 2.0;         // A
constexpr double kNumberDensity = 0.02; // particles / A^3, liquid-scale

Lattice cubicLatticeForDensity(std::size_t particleCount, double numberDensity) {
    const double volume = static_cast<double>(particleCount) / numberDensity;
    const double length = std::cbrt(volume);
    std::array<std::array<double, 3>, 3> m{};
    m[0] = {length, 0.0, 0.0};
    m[1] = {0.0, length, 0.0};
    m[2] = {0.0, 0.0, length};
    return Lattice(m);
}

ParticleData randomParticles(std::size_t count, const Lattice& lattice, unsigned seed) {
    std::mt19937 rng(seed);
    std::uniform_real_distribution<double> unit(0.0, 1.0);
    ParticleData particles;
    particles.resize(count);
    for (std::size_t i = 0; i < count; ++i) {
        const auto cart = lattice.fractionalToCartesian({unit(rng), unit(rng), unit(rng)});
        particles.x[i] = cart[0];
        particles.y[i] = cart[1];
        particles.z[i] = cart[2];
    }
    return particles;
}

void BM_CellListBuild(benchmark::State& state) {
    const auto n = static_cast<std::size_t>(state.range(0));
    const Lattice lattice = cubicLatticeForDensity(n, kNumberDensity);
    const ParticleData particles = randomParticles(n, lattice, 12345);

    for (auto _ : state) {
        CellList cellList;
        cellList.build(particles, lattice, kCutoff + kSkin);
        benchmark::DoNotOptimize(cellList);
    }
    state.SetComplexityN(static_cast<benchmark::ComplexityN>(n));
    state.SetItemsProcessed(static_cast<std::int64_t>(state.iterations()) * static_cast<std::int64_t>(n));
}
BENCHMARK(BM_CellListBuild)
    ->Arg(1'000)
    ->Arg(10'000)
    ->Arg(100'000)
    ->Arg(1'000'000)
    ->Unit(benchmark::kMillisecond);

void BM_VerletListBuild(benchmark::State& state) {
    const auto n = static_cast<std::size_t>(state.range(0));
    const Lattice lattice = cubicLatticeForDensity(n, kNumberDensity);
    const ParticleData particles = randomParticles(n, lattice, 54321);

    for (auto _ : state) {
        VerletList list(kSkin);
        list.build(particles, lattice, kCutoff);
        benchmark::DoNotOptimize(list);
    }
    state.SetComplexityN(static_cast<benchmark::ComplexityN>(n));
    state.SetItemsProcessed(static_cast<std::int64_t>(state.iterations()) * static_cast<std::int64_t>(n));
}
BENCHMARK(BM_VerletListBuild)
    ->Arg(1'000)
    ->Arg(10'000)
    ->Arg(100'000)
    ->Arg(1'000'000)
    ->Unit(benchmark::kMillisecond);

void BM_VerletListPairsTraversal(benchmark::State& state) {
    const auto n = static_cast<std::size_t>(state.range(0));
    const Lattice lattice = cubicLatticeForDensity(n, kNumberDensity);
    const ParticleData particles = randomParticles(n, lattice, 99999);

    VerletList list(kSkin);
    list.build(particles, lattice, kCutoff);

    std::size_t pairCount = 0;
    for (auto _ : state) {
        const auto& pairs = list.pairs();
        pairCount = pairs.size();
        benchmark::DoNotOptimize(pairCount);
        benchmark::DoNotOptimize(pairs.data());
    }
    state.SetComplexityN(static_cast<benchmark::ComplexityN>(n));
    state.counters["pairs_per_particle"] =
        static_cast<double>(pairCount) * 2.0 / static_cast<double>(n);
}
BENCHMARK(BM_VerletListPairsTraversal)
    ->Arg(1'000)
    ->Arg(10'000)
    ->Arg(100'000)
    ->Arg(1'000'000)
    ->Unit(benchmark::kMillisecond);

} // namespace
