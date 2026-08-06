#include "engines/geometry_analysis/periodic_graph.hpp"

#include <algorithm>
#include <numeric>
#include <optional>
#include <utility>

namespace aleator::engines {

namespace {

using Offset = std::array<int, 3>;

Offset add(const Offset& a, const Offset& b) {
    return {a[0] + b[0], a[1] + b[1], a[2] + b[2]};
}

Offset negate(const Offset& a) { return {-a[0], -a[1], -a[2]}; }

bool isZero(const Offset& a) { return a[0] == 0 && a[1] == 0 && a[2] == 0; }

/// Union-find where each node's `offsetToParent` is the lattice-vector
/// offset such that unwrappedPosition(parent) = unwrappedPosition(node) +
/// offsetToParent (i.e. "how far you travel, in lattice-translation units,
/// walking from this node to its parent"). find() path-compresses AND
/// accumulates offsets so that after a call, offsetToParent[x] is always
/// the total offset from x directly to the (possibly new) root returned.
class OffsetUnionFind {
public:
    explicit OffsetUnionFind(std::size_t n)
        : parent_(n), offsetToParent_(n, Offset{0, 0, 0}), percolates_(n, false) {
        std::iota(parent_.begin(), parent_.end(), std::size_t{0});
    }

    /// Returns {root, offset from x to root}.
    std::pair<std::size_t, Offset> find(std::size_t x) {
        if (parent_[x] == x) {
            return {x, Offset{0, 0, 0}};
        }
        const auto [root, parentToRoot] = find(parent_[x]);
        const Offset xToRoot = add(offsetToParent_[x], parentToRoot);
        parent_[x] = root;
        offsetToParent_[x] = xToRoot;
        return {root, xToRoot};
    }

    /// Unites the sets containing `u` and `v`, where traversing the edge
    /// from `u` to `v` corresponds to offset `edgeOffset` (unwrappedPosition
    /// (v) = unwrappedPosition(u) + edgeOffset, in the two nodes' own
    /// locally-consistent unwrapped frame — see periodic_graph.hpp).
    /// Returns the net offset of the cycle this edge closes if `u` and `v`
    /// were already connected (a real cycle in the graph, not a new union);
    /// returns nullopt if they were in different sets (a tree edge, no
    /// cycle information yet).
    std::optional<Offset> unite(std::size_t u, std::size_t v, const Offset& edgeOffset) {
        const auto [rootU, offsetUToRootU] = find(u);
        const auto [rootV, offsetVToRootV] = find(v);
        if (rootU == rootV) {
            // Two paths from u to the common root: directly via
            // offsetUToRootU, and via this edge (u -> v, offset
            // edgeOffset) then v -> root (offsetVToRootV). The difference
            // between them is the net offset of the cycle this edge closes.
            const Offset viaEdge = add(edgeOffset, offsetVToRootV);
            const Offset cycleOffset{viaEdge[0] - offsetUToRootU[0], viaEdge[1] - offsetUToRootU[1],
                                      viaEdge[2] - offsetUToRootU[2]};
            if (!isZero(cycleOffset)) {
                percolates_[rootU] = true;
            }
            return cycleOffset;
        }
        // Attach rootV under rootU (arbitrary choice; no rank/size
        // heuristic needed at this graph scale). Need offsetToParent_[rootV]
        // such that unwrappedPosition(rootU) = unwrappedPosition(rootV) +
        // offsetToParent_[rootV].
        // unwrappedPosition(v) = unwrappedPosition(u) + edgeOffset (edge
        // definition); unwrappedPosition(u) = unwrappedPosition(rootU) +
        // 0 - ... expressed via the offsets already recorded:
        // unwrappedPosition(rootU) = unwrappedPosition(u) - offsetUToRootU
        //   is WRONG sign convention -- re-derive directly below instead of
        //   trusting algebra by inspection; offsetToParent semantics are
        //   fixed by find()'s accumulation, so mirror that exactly:
        // find() guarantees unwrappedPosition(root) = unwrappedPosition(x)
        //   + offsetToParent_[x] for the x it was called on (that's what
        //   the accumulation computes). So:
        //   unwrappedPosition(rootU) = unwrappedPosition(u) + offsetUToRootU
        //   unwrappedPosition(rootV) = unwrappedPosition(v) + offsetVToRootV
        //   unwrappedPosition(v) = unwrappedPosition(u) + edgeOffset
        // => unwrappedPosition(rootU)
        //      = unwrappedPosition(v) - edgeOffset + offsetUToRootU
        //      = (unwrappedPosition(rootV) - offsetVToRootV) - edgeOffset + offsetUToRootU
        //    so offsetToParent_[rootV] (= unwrappedPosition(rootU) -
        //    unwrappedPosition(rootV)) = offsetUToRootU - edgeOffset - offsetVToRootV
        const Offset rootVToRootU = add(offsetUToRootU, negate(add(edgeOffset, offsetVToRootV)));
        parent_[rootV] = rootU;
        offsetToParent_[rootV] = rootVToRootU;
        percolates_[rootU] = percolates_[rootU] || percolates_[rootV];
        return std::nullopt;
    }

    bool percolates(std::size_t x) { return percolates_[find(x).first]; }

private:
    std::vector<std::size_t> parent_;
    std::vector<Offset> offsetToParent_;
    std::vector<bool> percolates_;
};

} // namespace

std::vector<bool> accessibleComponents(std::size_t nodeCount, const std::vector<PeriodicEdge>& edges,
                                        double thresholdRadius) {
    OffsetUnionFind uf(nodeCount);
    for (const PeriodicEdge& e : edges) {
        if (e.weight >= thresholdRadius) {
            uf.unite(e.from, e.to, e.imageOffset);
        }
    }
    std::vector<bool> result(nodeCount);
    for (std::size_t i = 0; i < nodeCount; ++i) {
        result[i] = uf.percolates(i);
    }
    return result;
}

double percolationThreshold(std::size_t nodeCount, const std::vector<PeriodicEdge>& edges) {
    std::vector<std::size_t> order(edges.size());
    std::iota(order.begin(), order.end(), std::size_t{0});
    std::sort(order.begin(), order.end(),
              [&](std::size_t a, std::size_t b) { return edges[a].weight > edges[b].weight; });

    OffsetUnionFind uf(nodeCount);
    for (std::size_t idx : order) {
        const PeriodicEdge& e = edges[idx];
        const auto cycleOffset = uf.unite(e.from, e.to, e.imageOffset);
        if (cycleOffset.has_value() && !isZero(*cycleOffset)) {
            return e.weight;
        }
    }
    return 0.0;
}

} // namespace aleator::engines
