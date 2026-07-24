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

/// Philox-4x32-10 (Salmon, Moraes, Dror & Shaw 2011), the reference
/// counter-based RNG for this codebase. Declared now so engines/ can be
/// written against CounterBasedRng; the Philox round function itself is a
/// correctness-critical piece of math and is deliberately not implemented
/// yet (CLAUDE.md invariant #1: no physics/math kernel ships before its
/// validation test). Calling any method throws NotImplemented.
class Philox4x32Rng final : public CounterBasedRng {
public:
    void seed(std::uint64_t seedValue, std::uint64_t streamIndex) override;
    std::uint64_t nextU64() override;
    double nextUniform() override;

private:
    // Not yet read or written by any implemented method (seed/nextU64/
    // nextUniform all throw NotImplemented) — [[maybe_unused]] silences
    // -Wunused-private-field until the Philox round function lands.
    [[maybe_unused]] std::array<std::uint32_t, 2> key_{};
    [[maybe_unused]] std::array<std::uint32_t, 4> counter_{};
};

} // namespace aleator::core
