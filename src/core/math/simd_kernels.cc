// Reference implementation of Highway runtime CPU dispatch. Structure
// mirrors Google's own hwy/examples/skeleton.cc exactly: this is the
// template every hot kernel in Aleator follows, so the pattern (rather than
// this particular kernel) is what matters here.

#include "core/math/simd_kernels.hpp"

// >>>> for dynamic dispatch only, skip if you want static dispatch

// First undef to prevent an error when this file is re-included by
// foreach_target.h itself.
#undef HWY_TARGET_INCLUDE
// Path is resolved against the include directories of this target (src/ is
// the include root — see cmake/AleatorTarget.cmake), matching how this file
// is referenced everywhere else: "core/math/simd_kernels.cc".
#define HWY_TARGET_INCLUDE "core/math/simd_kernels.cc"
#include "hwy/foreach_target.h" // IWYU pragma: keep

// <<<< end of dynamic dispatch

// Must come after foreach_target.h to avoid redefinition errors.
#include "hwy/highway.h"

HWY_BEFORE_NAMESPACE();
namespace aleator::core {
// Per-target namespace: required for dynamic dispatch so that the N_SSE4,
// N_AVX2, N_AVX3, ... instantiations of VectorSumImpl can co-exist in one
// translation unit without violating the one-definition rule.
namespace HWY_NAMESPACE {
namespace {

namespace hn = hwy::HWY_NAMESPACE;

std::int64_t VectorSumImpl(const std::int64_t* HWY_RESTRICT data, std::size_t count) {
    const hn::ScalableTag<std::int64_t> d;
    const std::size_t lanes = hn::Lanes(d);

    auto accumulator = hn::Zero(d);
    std::size_t i = 0;
    for (; i + lanes <= count; i += lanes) {
        accumulator = hn::Add(accumulator, hn::LoadU(d, data + i));
    }
    std::int64_t total = hn::ReduceSum(d, accumulator);

    for (; i < count; ++i) {
        total += data[i];
    }
    return total;
}

} // namespace
// NOLINTNEXTLINE(google-readability-namespace-comments)
} // namespace HWY_NAMESPACE
} // namespace aleator::core
HWY_AFTER_NAMESPACE();

// The dispatch table below must be compiled exactly once; foreach_target.h
// re-includes this whole file once per target, and HWY_ONCE is true for only
// the last of those passes.
#if HWY_ONCE

namespace aleator::core {

HWY_EXPORT(VectorSumImpl);

std::int64_t vectorSum(const std::int64_t* data, std::size_t count) {
    return HWY_DYNAMIC_DISPATCH(VectorSumImpl)(data, count);
}

} // namespace aleator::core

#endif // HWY_ONCE
