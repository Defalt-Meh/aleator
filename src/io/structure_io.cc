#include "io/structure_io.hpp"

#include "core/exceptions.hpp"

namespace aleator::io {

std::pair<core::ParticleData, core::Lattice> readCif(const std::filesystem::path&) {
    throw aleator::NotImplemented("readCif");
}

void writePdb(const std::filesystem::path&, const core::ParticleData&, const core::Lattice&) {
    throw aleator::NotImplemented("writePdb");
}

void writeLammpsData(const std::filesystem::path&, const core::ParticleData&,
                      const core::Lattice&) {
    throw aleator::NotImplemented("writeLammpsData");
}

void appendTrajectoryFrame(const std::filesystem::path&, const core::ParticleData&,
                            const core::Lattice&) {
    throw aleator::NotImplemented("appendTrajectoryFrame");
}

} // namespace aleator::io
