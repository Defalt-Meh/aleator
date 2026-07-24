#include "forcefield/parameters/parameter_loader.hpp"

#include "core/exceptions.hpp"

namespace aleator::forcefield {

ParameterSet loadParameterSet(const std::filesystem::path&) {
    throw aleator::NotImplemented("loadParameterSet");
}

} // namespace aleator::forcefield
