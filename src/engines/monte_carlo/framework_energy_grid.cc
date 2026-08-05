#include "engines/monte_carlo/framework_energy_grid.hpp"

#include <algorithm>
#include <cmath>
#include <stdexcept>

namespace aleator::engines {

namespace {
// Measured directly: an uncapped grid's mean interpolation error on the
// real IRMOF-1/CH4 system was ~1e17-1e18 K, entirely driven by nodes very
// close to a framework atom's core (raw 1/r^12 repulsion), even though the
// MEDIAN error was already small -- see
// tests/validation/test_framework_energy_grid.cc and
// FrameworkEnergyGrid::kEnergyCapKelvin's doc comment for why clamping is
// physically harmless, not a fudge.
double capped(double energy) { return std::min(energy, FrameworkEnergyGrid::kEnergyCapKelvin); }
} // namespace

FrameworkEnergyGrid::FrameworkEnergyGrid(const forcefield::LennardJones& forceField,
                                          const core::ParticleData& frameworkParticles,
                                          const core::Lattice& lattice,
                                          const std::vector<std::uint32_t>& guestSpeciesIds,
                                          double spacingAngstrom)
    : lattice_(lattice), spacingAngstrom_(spacingAngstrom) {
    if (!(spacingAngstrom_ > 0.0)) {
        throw std::invalid_argument("FrameworkEnergyGrid: spacingAngstrom must be positive");
    }
    // Sized off perpendicularWidth (Å per lattice direction), same quantity
    // core::CellList sizes its grid from — round rather than floor (unlike
    // CellList) since accuracy, not a candidate-superset guarantee, is what
    // matters here: rounding keeps the achieved spacing close to what was
    // asked for in either direction.
    for (int axis = 0; axis < 3; ++axis) {
        const double width = lattice_.perpendicularWidth(axis);
        const int n = static_cast<int>(std::lround(width / spacingAngstrom_));
        pointCounts_[static_cast<std::size_t>(axis)] = n < 1 ? 1 : n;
    }

    std::vector<std::uint32_t> uniqueSpecies = guestSpeciesIds;
    std::sort(uniqueSpecies.begin(), uniqueSpecies.end());
    uniqueSpecies.erase(std::unique(uniqueSpecies.begin(), uniqueSpecies.end()), uniqueSpecies.end());

    const auto totalPoints = static_cast<std::size_t>(pointCounts_[0]) *
                              static_cast<std::size_t>(pointCounts_[1]) *
                              static_cast<std::size_t>(pointCounts_[2]);

    for (std::uint32_t guestSpecies : uniqueSpecies) {
        Table table;
        table.values.resize(totalPoints);
        for (int iz = 0; iz < pointCounts_[2]; ++iz) {
            for (int iy = 0; iy < pointCounts_[1]; ++iy) {
                for (int ix = 0; ix < pointCounts_[0]; ++ix) {
                    const std::array<double, 3> fractional{
                        static_cast<double>(ix) / static_cast<double>(pointCounts_[0]),
                        static_cast<double>(iy) / static_cast<double>(pointCounts_[1]),
                        static_cast<double>(iz) / static_cast<double>(pointCounts_[2])};
                    const auto cart = lattice_.fractionalToCartesian(fractional);

                    double energy = 0.0;
                    for (std::size_t f = 0; f < frameworkParticles.size(); ++f) {
                        const auto params = forceField.mixedParameters(
                            guestSpecies, frameworkParticles.species[f]);
                        const auto d = lattice_.minimumImageDisplacement(
                            cart, {frameworkParticles.x[f], frameworkParticles.y[f],
                                   frameworkParticles.z[f]});
                        const double r = std::sqrt(d[0] * d[0] + d[1] * d[1] + d[2] * d[2]);
                        energy += forceField.pairEnergy(r, params);
                    }
                    table.values[flatten(ix, iy, iz)] = capped(energy);
                }
            }
        }
        tables_.emplace(guestSpecies, std::move(table));
    }
}

double FrameworkEnergyGrid::interpolate(std::uint32_t guestSpeciesId,
                                         const std::array<double, 3>& position) const {
    const Table& table = tables_.at(guestSpeciesId);

    auto fractional = lattice_.cartesianToFractional(position);
    for (double& f : fractional) {
        f -= std::floor(f);
    }

    std::array<int, 3> i0{};
    std::array<int, 3> i1{};
    std::array<double, 3> t{};
    for (int axis = 0; axis < 3; ++axis) {
        const int n = pointCounts_[static_cast<std::size_t>(axis)];
        const double g = fractional[static_cast<std::size_t>(axis)] * static_cast<double>(n);
        const int base = static_cast<int>(std::floor(g));
        t[static_cast<std::size_t>(axis)] = g - static_cast<double>(base);
        const int wrapped = ((base % n) + n) % n;
        i0[static_cast<std::size_t>(axis)] = wrapped;
        i1[static_cast<std::size_t>(axis)] = (wrapped + 1) % n;
    }

    auto value = [&](int ix, int iy, int iz) { return table.values[flatten(ix, iy, iz)]; };

    const double c00 =
        value(i0[0], i0[1], i0[2]) * (1.0 - t[0]) + value(i1[0], i0[1], i0[2]) * t[0];
    const double c10 =
        value(i0[0], i1[1], i0[2]) * (1.0 - t[0]) + value(i1[0], i1[1], i0[2]) * t[0];
    const double c01 =
        value(i0[0], i0[1], i1[2]) * (1.0 - t[0]) + value(i1[0], i0[1], i1[2]) * t[0];
    const double c11 =
        value(i0[0], i1[1], i1[2]) * (1.0 - t[0]) + value(i1[0], i1[1], i1[2]) * t[0];

    const double c0 = c00 * (1.0 - t[1]) + c10 * t[1];
    const double c1 = c01 * (1.0 - t[1]) + c11 * t[1];

    return c0 * (1.0 - t[2]) + c1 * t[2];
}

} // namespace aleator::engines
