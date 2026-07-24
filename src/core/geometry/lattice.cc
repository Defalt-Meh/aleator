#include "core/geometry/lattice.hpp"

#include "core/exceptions.hpp"

namespace aleator::core {

std::array<double, 3> Lattice::fractionalToCartesian(const std::array<double, 3>&) const {
    throw NotImplemented("Lattice::fractionalToCartesian");
}

std::array<double, 3> Lattice::cartesianToFractional(const std::array<double, 3>&) const {
    throw NotImplemented("Lattice::cartesianToFractional");
}

std::array<double, 3> Lattice::minimumImageDisplacement(const std::array<double, 3>&,
                                                          const std::array<double, 3>&) const {
    throw NotImplemented("Lattice::minimumImageDisplacement");
}

std::array<double, 3> Lattice::wrapIntoCell(const std::array<double, 3>&) const {
    throw NotImplemented("Lattice::wrapIntoCell");
}

} // namespace aleator::core
