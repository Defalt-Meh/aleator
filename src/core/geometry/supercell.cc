#include "core/geometry/supercell.hpp"

#include <stdexcept>

namespace aleator::core {

std::array<int, 3> minimumSupercellReplication(const Lattice& lattice, double cutoffAngstrom) {
    if (!(cutoffAngstrom > 0.0)) {
        throw std::invalid_argument("minimumSupercellReplication: cutoffAngstrom must be positive");
    }
    std::array<int, 3> n{1, 1, 1};
    for (int axis = 0; axis < 3; ++axis) {
        const double width = lattice.perpendicularWidth(axis);
        while (static_cast<double>(n[static_cast<std::size_t>(axis)]) * width <=
               2.0 * cutoffAngstrom) {
            ++n[static_cast<std::size_t>(axis)];
        }
    }
    return n;
}

SupercellResult replicateSupercell(const ParticleData& particles, const Lattice& lattice, int nx,
                                    int ny, int nz) {
    if (nx < 1 || ny < 1 || nz < 1) {
        throw std::invalid_argument("replicateSupercell: nx, ny, nz must each be >= 1");
    }

    const auto& m = lattice.matrix();
    const std::size_t n = particles.size();
    const auto totalCells = static_cast<std::size_t>(nx) * static_cast<std::size_t>(ny) *
                             static_cast<std::size_t>(nz);

    ParticleData replicated;
    replicated.reserve(n * totalCells);

    for (int i = 0; i < nx; ++i) {
        for (int j = 0; j < ny; ++j) {
            for (int k = 0; k < nz; ++k) {
                const double ox = i * m[0][0] + j * m[1][0] + k * m[2][0];
                const double oy = i * m[0][1] + j * m[1][1] + k * m[2][1];
                const double oz = i * m[0][2] + j * m[1][2] + k * m[2][2];
                for (std::size_t p = 0; p < n; ++p) {
                    replicated.push_back(particles.x[p] + ox, particles.y[p] + oy,
                                          particles.z[p] + oz, particles.vx[p], particles.vy[p],
                                          particles.vz[p], particles.mass[p], particles.charge[p],
                                          particles.species[p]);
                }
            }
        }
    }

    std::array<std::array<double, 3>, 3> newMatrix{};
    for (int d = 0; d < 3; ++d) {
        newMatrix[0][static_cast<std::size_t>(d)] = nx * m[0][static_cast<std::size_t>(d)];
        newMatrix[1][static_cast<std::size_t>(d)] = ny * m[1][static_cast<std::size_t>(d)];
        newMatrix[2][static_cast<std::size_t>(d)] = nz * m[2][static_cast<std::size_t>(d)];
    }

    return SupercellResult{std::move(replicated), Lattice(newMatrix)};
}

} // namespace aleator::core
