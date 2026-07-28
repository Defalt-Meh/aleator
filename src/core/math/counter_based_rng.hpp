#pragma once

#include <array>
#include <cstdint>

namespace aleator::core {

/// Interface for a counter-based RNG (CLAUDE.md invariant #5: seeded,
/// counter-based — Philox/Threefry style — one independent stream per
/// thread; same seed + same input + same thread count must give a
/// bit-identical trajectory).
///
/// Unlike a linear/LCG generator, a counter-based generator maps
/// (key, counter) -> pseudorandom output with no mutable internal state
/// beyond the counter, which is what makes it possible to give every thread
/// (or every MC trial move) its own independent, reproducible stream by
/// construction rather than by careful ordering of a single shared stream.
class CounterBasedRng {
public:
    virtual ~CounterBasedRng() = default;

    /// (seedValue, streamIndex) together select the key. Two RNGs
    /// constructed with the same pair produce identical output regardless
    /// of what else is running concurrently.
    virtual void seed(std::uint64_t seedValue, std::uint64_t streamIndex) = 0;

    /// Advances the counter and returns the next 64-bit output.
    virtual std::uint64_t nextU64() = 0;

    /// Returns a uniform double in [0, 1).
    virtual double nextUniform() = 0;
};

/// The raw Philox4x32-10 bijection (Salmon, Moraes, Dror & Shaw 2011,
/// "Parallel Random Numbers: As Easy as 1, 2, 3"): maps a 128-bit counter
/// and 64-bit key to a 128-bit output through 10 rounds of a Feistel-like
/// mix, with no higher-level "stream" semantics layered on top. Exposed
/// directly (rather than only through Philox4x32Rng) so it can be validated
/// against the official Random123 known-answer-test vectors bit-for-bit
/// (CLAUDE.md invariant #1: reference test before implementation) — see
/// tests/validation/test_philox_random123_vectors.cc and counter_based_rng.cc
/// for the algorithm itself.
[[nodiscard]] std::array<std::uint32_t, 4> philox4x32_10(std::array<std::uint32_t, 4> counter,
                                                           std::array<std::uint32_t, 2> key);

/// Philox4x32-10 (Salmon, Moraes, Dror & Shaw 2011), the reference
/// counter-based RNG for this codebase. seed(seedValue, streamIndex) sets
/// the 64-bit key to seedValue and the upper 64 bits of the counter to
/// streamIndex, permanently separating any two streams with a different
/// (seedValue, streamIndex) pair; the lower 64 bits of the counter are an
/// internal call index that advances by one on every nextU64(), which is
/// what gives repeated calls on one instance a deterministic, non-repeating
/// sequence. nextU64() uses only the first two of the four 32-bit words
/// philox4x32_10 produces per call — simpler and correctness-first, not
/// yet a target for optimization (CLAUDE.md #5: performance work happens
/// only after correctness, and only where the profiler points).
class Philox4x32Rng final : public CounterBasedRng {
public:
    void seed(std::uint64_t seedValue, std::uint64_t streamIndex) override;
    std::uint64_t nextU64() override;
    double nextUniform() override;

private:
    std::array<std::uint32_t, 2> key_{};
    std::array<std::uint32_t, 4> counter_{}; // counter_[2..3] = streamIndex, fixed for this instance
    std::uint64_t callIndex_ = 0;
};

} // namespace aleator::core
