// Verifies the Highway dynamic-dispatch template (core/math/simd_kernels.*)
// produces identical results on every target Highway compiled AND that this
// CPU actually supports. hwy::SetSupportedTargetsForTest is the
// Highway-documented way to force HWY_DYNAMIC_DISPATCH through each target
// in turn without needing gtest or hand-rolled per-namespace calls.

#include <catch2/catch_test_macros.hpp>

#include <cstdint>
#include <numeric>
#include <vector>

#include <hwy/highway.h>

#include "core/math/simd_kernels.hpp"

namespace {

std::int64_t scalarSum(const std::vector<std::int64_t>& data) {
    return std::accumulate(data.begin(), data.end(), std::int64_t{0});
}

} // namespace

TEST_CASE("vectorSum matches a scalar reference on every compiled and supported target",
          "[unit][simd]") {
    // Deliberately not a multiple of any lane count, to exercise the
    // scalar remainder loop on every target.
    std::vector<std::int64_t> data(137);
    std::iota(data.begin(), data.end(), -50);
    const std::int64_t expected = scalarSum(data);

    const int64_t supported = hwy::SupportedTargets();
    REQUIRE(supported != 0);

    for (int64_t remaining = supported; remaining != 0; remaining &= (remaining - 1)) {
        const int64_t oneTarget = remaining & (-remaining);
        hwy::SetSupportedTargetsForTest(oneTarget);
        INFO("target: " << hwy::TargetName(oneTarget));
        REQUIRE(aleator::core::vectorSum(data.data(), data.size()) == expected);
    }
    hwy::SetSupportedTargetsForTest(0); // restore normal CPU-detected behavior
}

TEST_CASE("vectorSum handles empty and single-element input", "[unit][simd]") {
    REQUIRE(aleator::core::vectorSum(nullptr, 0) == 0);

    const std::int64_t single = 42;
    REQUIRE(aleator::core::vectorSum(&single, 1) == 42);
}
