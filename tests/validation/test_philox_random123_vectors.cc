// Validation anchor (CLAUDE.md invariant #5: deterministic, reproducible,
// counter-based RNG). Philox4x32-10 is a fixed, published bijection, so
// unlike a validation against a physics simulation there is a hard,
// bit-exact reference: the official known-answer-test vectors from the
// Random123 reference implementation (DEShawResearch/random123,
// tests/kat_vectors, fetched 2026-07-28), reproduced verbatim below. If
// this doesn't match bit-for-bit, the implementation is wrong — there is no
// tolerance to loosen (CLAUDE.md invariant #1).

#include <catch2/catch_test_macros.hpp>

#include <array>
#include <cstdint>

#include "core/math/counter_based_rng.hpp"

using aleator::core::philox4x32_10;

TEST_CASE("philox4x32_10 matches the official Random123 known-answer-test vectors",
          "[validation][rng]") {
    // "philox4x32 10 <ctr0> <ctr1> <ctr2> <ctr3> <key0> <key1>  <out0> <out1> <out2> <out3>"
    SECTION("all-zero counter and key") {
        const auto out = philox4x32_10({0x00000000u, 0x00000000u, 0x00000000u, 0x00000000u},
                                        {0x00000000u, 0x00000000u});
        const std::array<std::uint32_t, 4> expected{0x6627e8d5u, 0xe169c58du, 0xbc57ac4cu,
                                                      0x9b00dbd8u};
        REQUIRE(out == expected);
    }

    SECTION("all-ones counter and key") {
        const auto out = philox4x32_10({0xffffffffu, 0xffffffffu, 0xffffffffu, 0xffffffffu},
                                        {0xffffffffu, 0xffffffffu});
        const std::array<std::uint32_t, 4> expected{0x408f276du, 0x41c83b0eu, 0xa20bc7c6u,
                                                      0x6d5451fdu};
        REQUIRE(out == expected);
    }

    SECTION("pi-digit counter and key") {
        const auto out = philox4x32_10({0x243f6a88u, 0x85a308d3u, 0x13198a2eu, 0x03707344u},
                                        {0xa4093822u, 0x299f31d0u});
        const std::array<std::uint32_t, 4> expected{0xd16cfe09u, 0x94fdcceb, 0x5001e420u,
                                                      0x24126ea1u};
        REQUIRE(out == expected);
    }
}

TEST_CASE("Philox4x32Rng is deterministic per (seed, streamIndex) and streams don't collide",
          "[validation][rng]") {
    aleator::core::Philox4x32Rng a;
    a.seed(42, 0);
    aleator::core::Philox4x32Rng b;
    b.seed(42, 0);

    for (int i = 0; i < 8; ++i) {
        REQUIRE(a.nextU64() == b.nextU64());
    }

    aleator::core::Philox4x32Rng streamZero;
    streamZero.seed(42, 0);
    aleator::core::Philox4x32Rng streamOne;
    streamOne.seed(42, 1);
    // Different streamIndex -> different key/counter entirely -> the first
    // draw must not coincide (astronomically unlikely by chance for a
    // 64-bit output; a collision here would indicate streamIndex is not
    // actually wired into the counter/key).
    REQUIRE(streamZero.nextU64() != streamOne.nextU64());

    aleator::core::Philox4x32Rng differentSeed;
    differentSeed.seed(43, 0);
    aleator::core::Philox4x32Rng originalAgain;
    originalAgain.seed(42, 0);
    REQUIRE(differentSeed.nextU64() != originalAgain.nextU64());
}

TEST_CASE("Philox4x32Rng::nextUniform stays within [0, 1) over many draws", "[validation][rng]") {
    aleator::core::Philox4x32Rng rng;
    rng.seed(7, 3);
    double sum = 0.0;
    constexpr int kDraws = 100000;
    for (int i = 0; i < kDraws; ++i) {
        const double u = rng.nextUniform();
        REQUIRE(u >= 0.0);
        REQUIRE(u < 1.0);
        sum += u;
    }
    // Mean of kDraws uniform(0,1) draws should be close to 0.5; standard
    // error of the mean here is 1/sqrt(12*kDraws) ~ 9.1e-4, so a 5-sigma
    // band is a generous, non-flaky bound on a genuine implementation while
    // still catching a badly-biased one (e.g. a byte-order or shift bug).
    const double mean = sum / kDraws;
    REQUIRE(mean > 0.5 - 5.0 * 9.13e-4);
    REQUIRE(mean < 0.5 + 5.0 * 9.13e-4);
}
