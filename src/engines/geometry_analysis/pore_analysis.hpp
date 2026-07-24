#pragma once

#include "core/geometry/lattice.hpp"
#include "core/math/particle_data.hpp"

namespace aleator::engines {

/// Porous-material geometry descriptors (Voronoi-based pore limiting
/// diameter / largest cavity diameter / accessible surface area), validated
/// against Zeo++'s published diameters for LTA, MFI, FAU (CLAUDE.md #4,
/// tolerance 1e-3 Å) before being trusted. `particles` are treated as rigid
/// framework atoms with van-der-Waals radii supplied by the caller via the
/// (not-yet-existing) force-field parameter set; this header only declares
/// the entry points.
struct PoreGeometry {
    double largestCavityDiameterAngstrom;
    double poreLimitingDiameterAngstrom;
    double accessibleSurfaceAreaAngstromSq;
};

/// Throws NotImplemented.
[[nodiscard]] PoreGeometry analyzePoreGeometry(const core::ParticleData& framework,
                                                const core::Lattice& lattice);

} // namespace aleator::engines
