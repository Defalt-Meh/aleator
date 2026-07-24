#include "core/math/counter_based_rng.hpp"

#include "core/exceptions.hpp"

namespace aleator::core {

void Philox4x32Rng::seed(std::uint64_t, std::uint64_t) {
    throw NotImplemented("Philox4x32Rng::seed");
}

std::uint64_t Philox4x32Rng::nextU64() { throw NotImplemented("Philox4x32Rng::nextU64"); }

double Philox4x32Rng::nextUniform() { throw NotImplemented("Philox4x32Rng::nextUniform"); }

} // namespace aleator::core
