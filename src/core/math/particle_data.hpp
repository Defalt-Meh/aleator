#pragma once

#include <cstddef>
#include <cstdint>
#include <vector>

#include "core/memory/aligned_allocator.hpp"

namespace aleator::core {

/// Structure-of-arrays particle storage (CLAUDE.md invariant #2: SoA,
/// always — never `std::vector<Particle>`).
///
/// This is a plain data container: resizing and indexing only. It carries no
/// physics (no energy/force evaluation, no dynamics) — those live in
/// forcefield/ and engines/ and operate on a ParticleData by reference.
///
/// Units (CLAUDE.md #2.4): position in Å, velocity in Å/ps, mass in amu,
/// charge in e.
struct ParticleData {
    std::vector<double, AlignedAllocator<double>> x, y, z;
    std::vector<double, AlignedAllocator<double>> vx, vy, vz;
    std::vector<double, AlignedAllocator<double>> mass;
    std::vector<double, AlignedAllocator<double>> charge;
    std::vector<std::uint32_t> species;

    [[nodiscard]] std::size_t size() const noexcept { return x.size(); }
    [[nodiscard]] bool empty() const noexcept { return x.empty(); }

    void reserve(std::size_t n) {
        x.reserve(n);
        y.reserve(n);
        z.reserve(n);
        vx.reserve(n);
        vy.reserve(n);
        vz.reserve(n);
        mass.reserve(n);
        charge.reserve(n);
        species.reserve(n);
    }

    void resize(std::size_t n) {
        x.resize(n);
        y.resize(n);
        z.resize(n);
        vx.resize(n);
        vy.resize(n);
        vz.resize(n);
        mass.resize(n);
        charge.resize(n);
        species.resize(n);
    }

    /// Appends one particle at the position/velocity/mass/charge/species
    /// given. Returns the index of the newly added particle.
    std::size_t push_back(double px, double py, double pz, double pvx, double pvy, double pvz,
                           double pmass, double pcharge, std::uint32_t pspecies) {
        x.push_back(px);
        y.push_back(py);
        z.push_back(pz);
        vx.push_back(pvx);
        vy.push_back(pvy);
        vz.push_back(pvz);
        mass.push_back(pmass);
        charge.push_back(pcharge);
        species.push_back(pspecies);
        return x.size() - 1;
    }

    /// Overwrites particle `dest`'s fields with particle `src`'s. Used by
    /// GCMC deletion (engines/monte_carlo) to implement O(1) removal:
    /// copy the last particle into the removed slot, then popBack().
    void copyParticle(std::size_t dest, std::size_t src) {
        x[dest] = x[src];
        y[dest] = y[src];
        z[dest] = z[src];
        vx[dest] = vx[src];
        vy[dest] = vy[src];
        vz[dest] = vz[src];
        mass[dest] = mass[src];
        charge[dest] = charge[src];
        species[dest] = species[src];
    }

    /// Drops the last particle. Paired with copyParticle() for O(1)
    /// swap-and-pop removal by index (order of remaining particles is not
    /// preserved).
    void popBack() {
        x.pop_back();
        y.pop_back();
        z.pop_back();
        vx.pop_back();
        vy.pop_back();
        vz.pop_back();
        mass.pop_back();
        charge.pop_back();
        species.pop_back();
    }
};

} // namespace aleator::core
