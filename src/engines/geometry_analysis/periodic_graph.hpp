#pragma once

#include <array>
#include <cstddef>
#include <vector>

namespace aleator::engines {

/// One edge of an abstract graph embedded in a periodic (crystallographic)
/// space: `from`/`to` are node indices, `weight` is the "radius" available
/// along this edge (the constraining quantity for percolation — e.g. a
/// Voronoi-network edge's free-sphere radius), and `imageOffset` is the net
/// number of lattice-vector translations accumulated when physically
/// traversing this edge from `from` to `to` (see PeriodicEdge's use in
/// voronoi_network.hpp for how this is derived from real geometry).
/// (0, 0, 0) means the edge stays within a single reference image; nonzero
/// means it crosses a periodic boundary.
struct PeriodicEdge {
    std::size_t from;
    std::size_t to;
    double weight;
    std::array<int, 3> imageOffset;
};

/// Largest threshold T such that the subgraph of edges with weight >= T
/// contains a cycle with nonzero net imageOffset — i.e. the largest probe
/// radius that can still percolate through the periodic structure (escape
/// its own reference image and return to an equivalent node by a genuinely
/// different path). Returns 0.0 if no such threshold exists for any
/// non-negative weight (the graph never percolates, e.g. a fully isolated
/// pocket).
///
/// Algorithm: process edges in descending weight order with a union-find
/// that additionally tracks each node's cumulative imageOffset relative to
/// its set's root (a standard "weighted/potentialized union-find," the same
/// technique used to detect a graph's cycle space parity in other domains).
/// Union two different sets: attach normally, recording the offset needed to
/// keep the invariant consistent. Encounter an edge between two nodes
/// ALREADY in the same set: this closes a cycle; the cycle's net offset is
/// computable from the two nodes' current root-relative offsets, and if that
/// net offset is nonzero, the graph percolates at this edge's weight (the
/// first — i.e. largest-weight — such edge in the descending sweep gives the
/// answer, since threshold subgraphs only gain edges as the threshold
/// decreases).
[[nodiscard]] double percolationThreshold(std::size_t nodeCount,
                                           const std::vector<PeriodicEdge>& edges);

/// For a fixed probe radius (rather than percolationThreshold()'s "largest
/// radius at which percolation exists anywhere"): partitions the subgraph
/// of edges with weight >= thresholdRadius into connected components (an
/// edge's weight >= thresholdRadius implies both endpoints' own radii do
/// too, by construction — see voronoi_network.hpp), and returns, per node,
/// whether that node's component percolates (a nonzero-offset cycle exists
/// within it) — the exact/graph-based classification behind
/// PoreGeometry's accessible-vs-inaccessible split (Zeo++'s "channel" vs
/// "pocket"). A node whose own radius is below thresholdRadius keeps
/// whatever component it starts in (a singleton, trivially non-percolating
/// — such nodes are never usable by a probe this size to begin with, so
/// their accessibility flag is never meaningfully queried by a caller).
[[nodiscard]] std::vector<bool> accessibleComponents(std::size_t nodeCount,
                                                      const std::vector<PeriodicEdge>& edges,
                                                      double thresholdRadius);

} // namespace aleator::engines
