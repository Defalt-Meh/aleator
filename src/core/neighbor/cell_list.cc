#include "core/neighbor/cell_list.hpp"

#include <array>
#include <cmath>
#include <stdexcept>

namespace aleator::core {

namespace {

/// Distinct raw offsets along one axis, after collapsing offsets that
/// wrap to the same target cell for a small cell count. For n>=3, {-1,0,1}
/// are pairwise distinct mod n and nothing collapses. For n==2, -1 and +1
/// both land on the same (only other) cell, so a single representative
/// (+1) stands in for both. For n==1, every offset lands on the only
/// cell, so only 0 is needed. Collapsing like this (rather than emitting
/// all of {-1,0,1} and deduplicating per-cell at scan time) is what keeps
/// forEachCandidatePair() from visiting the same source/target cell pair
/// more than once when the grid is smaller than 3 cells wide along some
/// axis — a real case (not just a theoretical edge case), since a cutoff
/// close to L_perp/2 forces exactly 2 cells along that axis.
std::array<int, 3> axisOffsets(int n, int& count) {
    if (n == 1) {
        count = 1;
        return {0, 0, 0};
    }
    if (n == 2) {
        count = 2;
        return {0, 1, 0};
    }
    count = 3;
    return {-1, 0, 1};
}

} // namespace

CellList::CellCoord CellList::wrap(CellCoord c, CellCoord offset) const noexcept {
    auto wrapAxis = [](int value, int n) {
        int m = value % n;
        return m < 0 ? m + n : m;
    };
    return {wrapAxis(c.x + offset.x, nx_), wrapAxis(c.y + offset.y, ny_), wrapAxis(c.z + offset.z, nz_)};
}

void CellList::build(const ParticleData& particles, const Lattice& lattice, double radius) {
    if (!(radius > 0.0)) {
        throw std::invalid_argument("CellList::build: radius must be positive");
    }

    auto cellCountForAxis = [&](int axis) {
        const int n = static_cast<int>(std::floor(lattice.perpendicularWidth(axis) / radius));
        return n < 1 ? 1 : n;
    };
    nx_ = cellCountForAxis(0);
    ny_ = cellCountForAxis(1);
    nz_ = cellCountForAxis(2);

    bins_.assign(static_cast<std::size_t>(nx_) * static_cast<std::size_t>(ny_) *
                     static_cast<std::size_t>(nz_),
                 {});

    auto clampBin = [](int bin, int n) {
        if (bin < 0) {
            return 0;
        }
        if (bin >= n) {
            return n - 1;
        }
        return bin;
    };

    for (std::size_t i = 0; i < particles.size(); ++i) {
        auto fractional =
            lattice.cartesianToFractional({particles.x[i], particles.y[i], particles.z[i]});
        for (double& f : fractional) {
            f -= std::floor(f);
        }
        const int bx = clampBin(static_cast<int>(fractional[0] * static_cast<double>(nx_)), nx_);
        const int by = clampBin(static_cast<int>(fractional[1] * static_cast<double>(ny_)), ny_);
        const int bz = clampBin(static_cast<int>(fractional[2] * static_cast<double>(nz_)), nz_);
        bins_[static_cast<std::size_t>(flatten({bx, by, bz}))].push_back(i);
    }

    int countX = 0;
    int countY = 0;
    int countZ = 0;
    const auto offsetsX = axisOffsets(nx_, countX);
    const auto offsetsY = axisOffsets(ny_, countY);
    const auto offsetsZ = axisOffsets(nz_, countZ);

    offsets_.clear();
    offsets_.reserve(static_cast<std::size_t>(countX * countY * countZ));
    for (int xi = 0; xi < countX; ++xi) {
        for (int yi = 0; yi < countY; ++yi) {
            for (int zi = 0; zi < countZ; ++zi) {
                offsets_.push_back({offsetsX[static_cast<std::size_t>(xi)],
                                     offsetsY[static_cast<std::size_t>(yi)],
                                     offsetsZ[static_cast<std::size_t>(zi)]});
            }
        }
    }
}

} // namespace aleator::core
