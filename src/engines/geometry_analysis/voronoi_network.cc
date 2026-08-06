#include "engines/geometry_analysis/voronoi_network.hpp"

#include <voro++.hh>

#include <algorithm>
#include <array>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <unordered_map>
#include <utility>

namespace aleator::engines {

namespace {

using Vec3 = std::array<double, 3>;

double dot(const Vec3& a, const Vec3& b) { return a[0] * b[0] + a[1] * b[1] + a[2] * b[2]; }
Vec3 sub(const Vec3& a, const Vec3& b) { return {a[0] - b[0], a[1] - b[1], a[2] - b[2]}; }
Vec3 add(const Vec3& a, const Vec3& b) { return {a[0] + b[0], a[1] + b[1], a[2] + b[2]}; }
Vec3 scale(const Vec3& a, double s) { return {a[0] * s, a[1] * s, a[2] * s}; }
double norm(const Vec3& a) { return std::sqrt(dot(a, a)); }

/// The six numbers Voro++'s periodic container requires: canonical
/// lower-triangular box vectors (bx,0,0), (bxy,by,0), (bxz,byz,bz), with
/// bx,by,bz > 0. Computed via Gram-Schmidt of `lattice`'s own vectors --
/// an isometric re-expression of the identical lattice (Gram-Schmidt
/// preserves the Gram/metric matrix, so cell lengths/angles/volume are
/// unchanged), not an approximation. A `core::Lattice` built from these six
/// numbers therefore has fractional coordinates that mean exactly the same
/// physical point, at exactly the same physical distances, as `lattice`'s
/// own fractional coordinates -- which is what lets buildVoronoiNetwork()
/// convert positions in and out by round-tripping through fractional
/// coordinates instead of carrying an explicit rotation matrix through the
/// whole computation.
struct CanonicalBox {
    double bx, bxy, by, bxz, byz, bz;
};

CanonicalBox computeCanonicalBox(const core::Lattice& lattice) {
    const auto& m = lattice.matrix();
    const Vec3 a0{m[0][0], m[0][1], m[0][2]};
    const Vec3 a1{m[1][0], m[1][1], m[1][2]};
    const Vec3 a2{m[2][0], m[2][1], m[2][2]};

    const double bx = norm(a0);
    const Vec3 r0 = scale(a0, 1.0 / bx);
    const double bxy = dot(a1, r0);
    const Vec3 a1Perp = sub(a1, scale(r0, bxy));
    const double by = norm(a1Perp);
    const Vec3 r1 = scale(a1Perp, 1.0 / by);
    const double bxz = dot(a2, r0);
    const double byz = dot(a2, r1);
    const Vec3 a2Perp = sub(a2, add(scale(r0, bxz), scale(r1, byz)));
    const double bz = norm(a2Perp);

    return {bx, bxy, by, bxz, byz, bz};
}

core::Lattice canonicalLatticeFrom(const CanonicalBox& box) {
    return core::Lattice(
        {{{box.bx, 0.0, 0.0}, {box.bxy, box.by, 0.0}, {box.bxz, box.byz, box.bz}}});
}

/// Spatial-hash key for deduplicating Voronoi vertices that different
/// cells report as (up to floating-point noise) the same geometric point.
/// 1e-6 Ang is many orders of magnitude above double-precision round-off at
/// atomic-structure scale and many orders below any real distinct-vertex
/// separation, so this cannot conflate genuinely different vertices.
constexpr double kQuantum = 1.0e-6;

struct QuantizedKey {
    long long x, y, z;
    bool operator==(const QuantizedKey& other) const noexcept {
        return x == other.x && y == other.y && z == other.z;
    }
};

struct QuantizedKeyHash {
    std::size_t operator()(const QuantizedKey& k) const noexcept {
        std::size_t h = std::hash<long long>()(k.x);
        h ^= std::hash<long long>()(k.y) + 0x9e3779b97f4a7c15ULL + (h << 6) + (h >> 2);
        h ^= std::hash<long long>()(k.z) + 0x9e3779b97f4a7c15ULL + (h << 6) + (h >> 2);
        return h;
    }
};

QuantizedKey keyFor(const Vec3& p) {
    return QuantizedKey{std::llround(p[0] / kQuantum), std::llround(p[1] / kQuantum),
                         std::llround(p[2] / kQuantum)};
}

} // namespace

VoronoiNetwork buildVoronoiNetwork(const core::ParticleData& framework, const core::Lattice& lattice,
                                    const std::vector<double>& perAtomRadiusAngstrom) {
    if (framework.size() != perAtomRadiusAngstrom.size()) {
        throw std::invalid_argument(
            "buildVoronoiNetwork: perAtomRadiusAngstrom size must match framework particle count");
    }
    if (framework.empty()) {
        throw std::invalid_argument("buildVoronoiNetwork: framework must not be empty");
    }

    const CanonicalBox box = computeCanonicalBox(lattice);
    const core::Lattice canonical = canonicalLatticeFrom(box);
    const std::size_t n = framework.size();

    // Voro++'s internal spatial-hash grid resolution: aim for roughly a
    // handful of particles per grid block, matching Voro++'s own examples'
    // sizing heuristic -- a performance tuning knob, not an accuracy one
    // (compute_cell() is exact regardless of grid resolution).
    const int perAxis = std::max(1, static_cast<int>(std::cbrt(static_cast<double>(n) / 4.0)));

    voro::container_periodic_poly con(box.bx, box.bxy, box.by, box.bxz, box.byz, box.bz, perAxis,
                                       perAxis, perAxis, 8);
    for (std::size_t i = 0; i < n; ++i) {
        auto frac = lattice.cartesianToFractional({framework.x[i], framework.y[i], framework.z[i]});
        for (double& f : frac) {
            f -= std::floor(f);
        }
        const auto canonicalPos = canonical.fractionalToCartesian(frac);
        con.put(static_cast<int>(i), canonicalPos[0], canonicalPos[1], canonicalPos[2],
                perAtomRadiusAngstrom[i]);
    }

    // Maps a Voro++ canonical-frame vertex position to (its physical
    // position back in `lattice`'s own frame, wrapped into the primary
    // cell) and (the integer lattice-translation this occurrence had been
    // offset by before wrapping) -- the latter is exactly a PeriodicEdge
    // imageOffset component once differenced between two vertices of the
    // same cell (see the edge-construction loop below).
    auto canonicalToOriginal = [&](const Vec3& canonicalPos) {
        auto frac = canonical.cartesianToFractional(
            {canonicalPos[0], canonicalPos[1], canonicalPos[2]});
        std::array<int, 3> imageIndex{};
        for (int d = 0; d < 3; ++d) {
            const double flo = std::floor(frac[static_cast<std::size_t>(d)]);
            imageIndex[static_cast<std::size_t>(d)] = static_cast<int>(flo);
            frac[static_cast<std::size_t>(d)] -= flo;
        }
        const auto originalPos = lattice.fractionalToCartesian(frac);
        return std::pair<Vec3, std::array<int, 3>>{Vec3{originalPos[0], originalPos[1], originalPos[2]},
                                                     imageIndex};
    };

    auto freeRadiusAt = [&](const Vec3& p) {
        double best = std::numeric_limits<double>::infinity();
        for (std::size_t j = 0; j < n; ++j) {
            const auto d =
                lattice.minimumImageDisplacement(p, {framework.x[j], framework.y[j], framework.z[j]});
            const double r =
                std::sqrt(d[0] * d[0] + d[1] * d[1] + d[2] * d[2]) - perAtomRadiusAngstrom[j];
            if (r < best) {
                best = r;
            }
        }
        return best;
    };

    VoronoiNetwork network;
    std::unordered_map<QuantizedKey, std::size_t, QuantizedKeyHash> nodeIndex;

    voro::c_loop_all_periodic vl(con);
    voro::voronoicell_neighbor cell;
    if (vl.start()) {
        do {
            double px = 0.0;
            double py = 0.0;
            double pz = 0.0;
            vl.pos(px, py, pz);
            if (!con.compute_cell(cell, vl)) {
                continue;
            }

            std::vector<double> verts;
            cell.vertices(px, py, pz, verts);
            const std::size_t vertexCount = verts.size() / 3;

            // Per-vertex-occurrence bookkeeping, local to this cell.
            std::vector<std::size_t> localNode(vertexCount);
            std::vector<std::array<int, 3>> localImageOffset(vertexCount);
            for (std::size_t v = 0; v < vertexCount; ++v) {
                const Vec3 raw{verts[3 * v], verts[3 * v + 1], verts[3 * v + 2]};
                auto [originalPos, imageIndex] = canonicalToOriginal(raw);
                const QuantizedKey key = keyFor(originalPos);
                const auto it = nodeIndex.find(key);
                std::size_t idx;
                if (it == nodeIndex.end()) {
                    idx = network.nodes.size();
                    network.nodes.push_back(
                        VoronoiNode{{originalPos[0], originalPos[1], originalPos[2]},
                                    freeRadiusAt(originalPos)});
                    nodeIndex.emplace(key, idx);
                } else {
                    idx = it->second;
                }
                localNode[v] = idx;
                localImageOffset[v] = imageIndex;
            }

            std::vector<int> faceVertices;
            cell.face_vertices(faceVertices);
            std::size_t pos = 0;
            while (pos < faceVertices.size()) {
                const auto count = static_cast<std::size_t>(faceVertices[pos]);
                ++pos;
                for (std::size_t k = 0; k < count; ++k) {
                    const auto vi = static_cast<std::size_t>(faceVertices[pos + k]);
                    const auto vj = static_cast<std::size_t>(faceVertices[pos + (k + 1) % count]);
                    const auto& offI = localImageOffset[vi];
                    const auto& offJ = localImageOffset[vj];
                    const std::array<int, 3> edgeOffset{offJ[0] - offI[0], offJ[1] - offI[1],
                                                          offJ[2] - offI[2]};

                    // Sample the free radius along the edge's actual raw
                    // (canonical-frame) segment, not just at its two
                    // endpoints — see VoronoiEdge's doc comment for why
                    // endpoints alone are provably insufficient (a real
                    // failure mode found while validating this milestone
                    // on LTA). Sampling density scales with edge length
                    // since edges range from ~1 Ang to (rarely, for atoms
                    // on special symmetry positions) most of a unit cell.
                    // kEdgeSampleSpacingAngstrom is chosen, not assumed: at
                    // 0.2 Ang, PLD already matched real Zeo++ output to
                    // <=0.003 Ang on LTA/MFI/FAU; at this 0.02 Ang (10x
                    // finer), it matches to <=0.0005 Ang on all three,
                    // comfortably inside CLAUDE.md's 1e-3 Ang validation
                    // target (see tests/validation/test_pore_geometry.cc for
                    // the exact numbers this was measured against). Real
                    // cost, not free: this is why pore analysis is a
                    // deliberately "run once per structure" computation
                    // (matching Zeo++'s own usage pattern), not a hot-loop
                    // operation.
                    constexpr double kEdgeSampleSpacingAngstrom = 0.02;
                    const Vec3 rawI{verts[3 * vi], verts[3 * vi + 1], verts[3 * vi + 2]};
                    const Vec3 rawJ{verts[3 * vj], verts[3 * vj + 1], verts[3 * vj + 2]};
                    const double edgeLength = norm(sub(rawJ, rawI));
                    const auto sampleCount = static_cast<std::size_t>(
                        std::max(2.0, std::ceil(edgeLength / kEdgeSampleSpacingAngstrom) + 1.0));
                    double minRadius = std::numeric_limits<double>::infinity();
                    for (std::size_t s = 0; s < sampleCount; ++s) {
                        const double t = static_cast<double>(s) / static_cast<double>(sampleCount - 1);
                        const Vec3 rawSample = add(rawI, scale(sub(rawJ, rawI), t));
                        const auto [samplePos, sampleImage] = canonicalToOriginal(rawSample);
                        (void)sampleImage;
                        const double r = freeRadiusAt(samplePos);
                        minRadius = std::min(minRadius, r);
                        network.maxFreeRadiusAlongAnyEdgeAngstrom =
                            std::max(network.maxFreeRadiusAlongAnyEdgeAngstrom, r);
                    }

                    network.edges.push_back(
                        VoronoiEdge{localNode[vi], localNode[vj], minRadius, edgeOffset});
                }
                pos += count;
            }
        } while (vl.inc());
    }

    return network;
}

} // namespace aleator::engines
