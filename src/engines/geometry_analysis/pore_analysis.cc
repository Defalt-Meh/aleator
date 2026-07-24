#include "engines/geometry_analysis/pore_analysis.hpp"

#include "core/exceptions.hpp"

namespace aleator::engines {

PoreGeometry analyzePoreGeometry(const core::ParticleData&, const core::Lattice&) {
    throw aleator::NotImplemented("analyzePoreGeometry");
}

} // namespace aleator::engines
