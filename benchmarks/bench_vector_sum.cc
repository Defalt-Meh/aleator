// Reference benchmark, paired with tests/unit/test_simd_kernels.cc. Every
// hot kernel gets one of these once it exists; this one exists purely to
// prove the harness (target, flags, JSON output) works end to end.

#include <cstdint>
#include <numeric>
#include <vector>

#include <benchmark/benchmark.h>

#include "core/math/simd_kernels.hpp"

namespace {

void BM_VectorSum(benchmark::State& state) {
    const auto count = static_cast<std::size_t>(state.range(0));
    std::vector<std::int64_t> data(count);
    std::iota(data.begin(), data.end(), 0);

    for (auto _ : state) {
        benchmark::DoNotOptimize(aleator::core::vectorSum(data.data(), data.size()));
    }
    const auto iterations = static_cast<std::int64_t>(state.iterations());
    const auto itemCount = static_cast<std::int64_t>(count);
    state.SetItemsProcessed(iterations * itemCount);
    state.SetBytesProcessed(iterations * itemCount * static_cast<std::int64_t>(sizeof(std::int64_t)));
}

} // namespace

BENCHMARK(BM_VectorSum)->Range(1 << 8, 1 << 20);
