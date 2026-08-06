// Unit tests for the abstract percolation-threshold algorithm
// (engines/geometry_analysis/periodic_graph.hpp), in isolation from any
// real Voronoi geometry -- synthetic graphs with known-by-construction
// percolation behavior, so the offset-tracking union-find's sign
// conventions are pinned down before it's trusted to compute a real pore
// limiting diameter.

#include <catch2/catch_test_macros.hpp>

#include <vector>

#include "engines/geometry_analysis/periodic_graph.hpp"

using aleator::engines::PeriodicEdge;
using aleator::engines::percolationThreshold;

TEST_CASE("percolationThreshold is 0 for a graph with no edges", "[unit][geometry]") {
    REQUIRE(percolationThreshold(3, {}) == 0.0);
}

TEST_CASE("percolationThreshold is 0 for a tree (no cycles at all)", "[unit][geometry]") {
    // 0 -- 1 -- 2, a path: no cycle exists regardless of offsets.
    const std::vector<PeriodicEdge> edges{
        {0, 1, 5.0, {0, 0, 0}},
        {1, 2, 3.0, {1, 0, 0}},
    };
    REQUIRE(percolationThreshold(3, edges) == 0.0);
}

TEST_CASE("percolationThreshold is 0 for a purely local cycle (zero net offset)",
          "[unit][geometry]") {
    // Triangle 0-1-2-0, all offsets zero: a real cycle, but it never leaves
    // the reference image, so it must not be mistaken for percolation.
    const std::vector<PeriodicEdge> edges{
        {0, 1, 5.0, {0, 0, 0}},
        {1, 2, 4.0, {0, 0, 0}},
        {2, 0, 3.0, {0, 0, 0}},
    };
    REQUIRE(percolationThreshold(3, edges) == 0.0);
}

TEST_CASE("percolationThreshold detects a single self-loop with nonzero offset",
          "[unit][geometry]") {
    // One node, one edge from it back to itself with a nonzero periodic
    // offset -- the degenerate but real case of a single atom per unit
    // cell in a simple lattice (verified directly against real Voro++
    // output during this milestone's design: a single-atom simple-cubic
    // cell's Voronoi faces are all self-periodic-neighbors).
    const std::vector<PeriodicEdge> edges{
        {0, 0, 5.0, {1, 0, 0}},
    };
    REQUIRE(percolationThreshold(1, edges) == 5.0);
}

TEST_CASE("percolationThreshold picks the bottleneck (minimum) weight along the "
          "percolating cycle, not the first or largest edge",
          "[unit][geometry]") {
    // 0 and 1 connected two ways: directly (offset 0, weight 3.0 -- a tree
    // edge) and via a periodic image (offset (1,0,0), weight 2.0). The
    // cycle these two edges form has net offset (1,0,0) -- nonzero, so it
    // percolates -- at the weaker (2.0) of the two edges.
    const std::vector<PeriodicEdge> edges{
        {0, 1, 3.0, {0, 0, 0}},
        {1, 0, 2.0, {1, 0, 0}},
    };
    REQUIRE(percolationThreshold(2, edges) == 2.0);
}

TEST_CASE("percolationThreshold ignores a weaker percolating cycle when a stronger "
          "one exists",
          "[unit][geometry]") {
    // Two independent percolating loops on disjoint node pairs: {0,1} with
    // bottleneck 2.0, {2,3} with bottleneck 4.0. The answer must be the
    // larger bottleneck (4.0) -- the whole point of processing edges in
    // descending weight order and stopping at the first cycle found.
    const std::vector<PeriodicEdge> edges{
        {0, 1, 3.0, {0, 0, 0}},
        {1, 0, 2.0, {1, 0, 0}},
        {2, 3, 5.0, {0, 0, 0}},
        {3, 2, 4.0, {0, 1, 0}},
    };
    REQUIRE(percolationThreshold(4, edges) == 4.0);
}

TEST_CASE("percolationThreshold requires a genuinely nonzero net offset, not just "
          "any second connection between components",
          "[unit][geometry]") {
    // 0 and 1 connected twice, both times with offset (0,0,0): this is
    // still just a local (non-percolating) cycle, however many parallel
    // edges connect the same pair.
    const std::vector<PeriodicEdge> edges{
        {0, 1, 5.0, {0, 0, 0}},
        {0, 1, 4.0, {0, 0, 0}},
    };
    REQUIRE(percolationThreshold(2, edges) == 0.0);
}

TEST_CASE("percolationThreshold detects percolation through a longer chain with a "
          "distant periodic closure",
          "[unit][geometry]") {
    // A chain 0-1-2-3-4 (all offset 0, i.e. these five nodes are all
    // "the same unwrapped position" up to the chain), closed by a single
    // edge 4->0 with a nonzero offset. The whole chain + closing edge
    // forms one percolating cycle; its bottleneck is the weakest edge in
    // it (1.5, the closing edge).
    const std::vector<PeriodicEdge> edges{
        {0, 1, 9.0, {0, 0, 0}}, {1, 2, 8.0, {0, 0, 0}}, {2, 3, 7.0, {0, 0, 0}},
        {3, 4, 6.0, {0, 0, 0}}, {4, 0, 1.5, {0, 0, 1}},
    };
    REQUIRE(percolationThreshold(5, edges) == 1.5);
}
