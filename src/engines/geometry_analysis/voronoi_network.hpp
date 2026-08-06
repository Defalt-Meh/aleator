#pragma once

#include <array>
#include <cstddef>
#include <vector>

#include "core/geometry/lattice.hpp"
#include "core/math/particle_data.hpp"

namespace aleator::engines {

/// One node of the periodic Voronoi network: a vertex of the radical
/// (radius-weighted) Voronoi decomposition of the framework, i.e. a
/// candidate center for a locally-largest empty sphere.
/// `freeRadiusAngstrom` is the TRUE free-sphere radius at this point (the
/// minimum, over every framework atom under minimum image, of
/// |position - atom| - atom.radius) — computed directly, not assumed equal
/// to any Voro++-internal quantity. This matters: Voro++ computes a power
/// (Laguerre) diagram, whose vertices are equidistant in *power distance*
/// (|p-c|^2 - r^2) to their generating atoms, which is not quite the same
/// locus as being equidistant in *free radius* (|p-c| - r) when those
/// atoms have different radii (the geometrically "exact" formulation for
/// polydisperse spheres is the Apollonius/additively-weighted diagram,
/// which is not planar-faced and not what any fast library computes). The
/// power diagram's vertices are nonetheless the standard, well-established
/// practical choice of candidate points (this is what Zeo++ itself does);
/// recomputing the true free radius at each candidate directly, rather
/// than trusting the power-diagram value, keeps this exact regardless of
/// that approximation's error.
struct VoronoiNode {
    std::array<double, 3> positionAngstrom; // wrapped into the primary cell
    double freeRadiusAngstrom;
};

/// One edge of the network: `from`/`to` are indices into `nodes`, connected
/// by a shared boundary of Voronoi cells (a probe can move continuously
/// between them). `minRadiusAlongEdgeAngstrom` is the minimum free-sphere
/// radius sampled along the edge's actual physical path — NOT just
/// min(endpoint radii): a real, measured failure mode found while
/// validating this milestone against real Zeo++ output is that some
/// Voronoi edges are anomalously long (observed: up to ~11.6 Å, spanning
/// almost an entire LTA unit cell, for an atom sitting on a special
/// symmetry position with a very open cell on one side) and pass directly
/// through a narrow window between two wide cavities — endpoint radii
/// alone completely miss that constriction, since both endpoints sit in
/// the wide regions on either side of it. Sampling along the edge's
/// length is what actually finds it (see voronoi_network.cc). `imageOffset`
/// is the number of lattice-vector translations accumulated when
/// physically traversing this edge from `from` to `to` (see
/// periodic_graph.hpp — this is exactly a PeriodicEdge's imageOffset).
struct VoronoiEdge {
    std::size_t from;
    std::size_t to;
    double minRadiusAlongEdgeAngstrom;
    std::array<int, 3> imageOffset;
};

struct VoronoiNetwork {
    std::vector<VoronoiNode> nodes;
    std::vector<VoronoiEdge> edges;
    /// The largest free-sphere radius found anywhere while sampling edge
    /// interiors (see VoronoiEdge) — not just at nodes/vertices. Tracked
    /// because it's a real, free (already-computed) improvement to the
    /// largest-cavity-diameter candidate set: the true maximal empty
    /// sphere for atoms of *different* radii is, strictly, a vertex of the
    /// additively-weighted (Apollonius) diagram, not the power/Laguerre
    /// diagram Voro++ actually computes, and those coincide only when the
    /// generating atoms have equal radii — so the true maximum can fall
    /// slightly off a power-diagram vertex, in a spot an edge sample can
    /// catch that a vertex-only search would miss. See
    /// tests/validation/test_pore_geometry.cc for how much this closes the
    /// gap against real Zeo++ output.
    double maxFreeRadiusAlongAnyEdgeAngstrom = 0.0;
};

/// Builds the periodic radical Voronoi network of `framework` under
/// `lattice`, via Voro++ (github.com/chr1shr/voro; CLAUDE.md section 0 /
/// README document why this depends on Voro++ rather than reimplementing
/// periodic weighted Voronoi decomposition). Exact for triclinic cells
/// (invariant #6): Voro++'s periodic container requires a specific
/// canonical lower-triangular basis, so this internally rotates into that
/// basis and rotates results back — a rotation (or, for a left-handed
/// input lattice, a reflection) is an isometry, so this changes nothing
/// about the physical distances/radii computed.
///
/// `perAtomRadiusAngstrom[i]` is the radius of `framework` particle i
/// (parallel array, same size as `framework`). Throws std::invalid_argument
/// if the sizes mismatch or if `framework` is empty.
[[nodiscard]] VoronoiNetwork buildVoronoiNetwork(
    const core::ParticleData& framework, const core::Lattice& lattice,
    const std::vector<double>& perAtomRadiusAngstrom);

} // namespace aleator::engines
