#pragma once

#include <cstddef>
#include <cstdint>

namespace aleator::core {

/// Sums `count` 64-bit integers starting at `data`. This is the reference
/// kernel for Highway runtime CPU dispatch (see simd_kernels.cc): every hot
/// kernel in this codebase follows this same HWY_EXPORT /
/// HWY_DYNAMIC_DISPATCH shape, compiled once per attainable target
/// (SSE4/AVX2/AVX-512/NEON/...) and selected at runtime by CPU detection.
/// No compile-time ISA assumption, no `-march=native` (CLAUDE.md #2.3, #6).
///
/// Integer (rather than floating-point) so that the result is exactly
/// reduction-order-independent: every compiled target, regardless of lane
/// width, must return the bit-identical sum (see
/// tests/unit/test_simd_kernels.cc). A floating-point kernel would not have
/// this property — see CLAUDE.md #10 on reassociation breaking
/// reproducibility.
///
/// `data` may be nullptr only if `count == 0`.
std::int64_t vectorSum(const std::int64_t* data, std::size_t count);

} // namespace aleator::core
