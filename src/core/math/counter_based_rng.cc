#include "core/math/counter_based_rng.hpp"

namespace aleator::core {

namespace {

// Philox4x32-10 constants, verbatim from the reference implementation
// (DEShawResearch/random123, include/Random123/philox.h): the two
// multipliers (chosen for good multiplicative-congruential statistical
// properties in the round function) and the two Weyl-sequence increments
// (2^32/golden-ratio and 2^32/sqrt(3), the standard choice ensuring the
// per-round key schedule never repeats within any realistic counter range).
constexpr std::uint32_t kPhiloxMul0 = 0xD2511F53u;
constexpr std::uint32_t kPhiloxMul1 = 0xCD9E8D57u;
constexpr std::uint32_t kPhiloxWeyl0 = 0x9E3779B9u;
constexpr std::uint32_t kPhiloxWeyl1 = 0xBB67AE85u;
constexpr int kPhiloxRounds = 10;

/// Splits the 64-bit product a*b into (hi, lo) 32-bit halves.
void mulhilo32(std::uint32_t a, std::uint32_t b, std::uint32_t& hi, std::uint32_t& lo) {
    const std::uint64_t product = static_cast<std::uint64_t>(a) * static_cast<std::uint64_t>(b);
    hi = static_cast<std::uint32_t>(product >> 32);
    lo = static_cast<std::uint32_t>(product);
}

/// One Philox4x32 round: ctr = (v0,v1,v2,v3), key = (k0,k1) ->
///   (hi1^v1^k0, lo1, hi0^v3^k1, lo0)
/// where (hi0,lo0) = mulhilo32(Mul0, v0), (hi1,lo1) = mulhilo32(Mul1, v2).
std::array<std::uint32_t, 4> philoxRound(const std::array<std::uint32_t, 4>& ctr,
                                          const std::array<std::uint32_t, 2>& key) {
    std::uint32_t hi0 = 0;
    std::uint32_t lo0 = 0;
    std::uint32_t hi1 = 0;
    std::uint32_t lo1 = 0;
    mulhilo32(kPhiloxMul0, ctr[0], hi0, lo0);
    mulhilo32(kPhiloxMul1, ctr[2], hi1, lo1);
    return {hi1 ^ ctr[1] ^ key[0], lo1, hi0 ^ ctr[3] ^ key[1], lo0};
}

} // namespace

std::array<std::uint32_t, 4> philox4x32_10(std::array<std::uint32_t, 4> counter,
                                            std::array<std::uint32_t, 2> key) {
    for (int round = 0; round < kPhiloxRounds; ++round) {
        if (round > 0) {
            // Key schedule: bump by the Weyl constants before every round
            // after the first (round 1 uses the original key).
            key[0] += kPhiloxWeyl0;
            key[1] += kPhiloxWeyl1;
        }
        counter = philoxRound(counter, key);
    }
    return counter;
}

void Philox4x32Rng::seed(std::uint64_t seedValue, std::uint64_t streamIndex) {
    key_ = {static_cast<std::uint32_t>(seedValue), static_cast<std::uint32_t>(seedValue >> 32)};
    counter_[0] = 0;
    counter_[1] = 0;
    counter_[2] = static_cast<std::uint32_t>(streamIndex);
    counter_[3] = static_cast<std::uint32_t>(streamIndex >> 32);
    callIndex_ = 0;
}

std::uint64_t Philox4x32Rng::nextU64() {
    counter_[0] = static_cast<std::uint32_t>(callIndex_);
    counter_[1] = static_cast<std::uint32_t>(callIndex_ >> 32);
    ++callIndex_;
    const auto out = philox4x32_10(counter_, key_);
    return (static_cast<std::uint64_t>(out[0]) << 32) | static_cast<std::uint64_t>(out[1]);
}

double Philox4x32Rng::nextUniform() {
    // Standard 53-bit-precision technique: top 53 bits of a 64-bit draw,
    // scaled into [0, 1). Not something invented for this codebase — the
    // same construction std::mt19937_64-based real distributions use.
    constexpr double kTwoPow53 = 9007199254740992.0; // 2^53
    return static_cast<double>(nextU64() >> 11) * (1.0 / kTwoPow53);
}

} // namespace aleator::core
