#include "engines/geometry_analysis/pore_analysis.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <limits>
#include <map>
#include <numbers>
#include <stdexcept>

#include "core/math/counter_based_rng.hpp"
#include "engines/geometry_analysis/periodic_graph.hpp"
#include "engines/geometry_analysis/voronoi_network.hpp"

namespace aleator::engines {

namespace {

using Vec3 = std::array<double, 3>;

double norm(const Vec3& v) { return std::sqrt(v[0] * v[0] + v[1] * v[1] + v[2] * v[2]); }

/// Every framework atom's radius, expanded from the per-species table.
std::vector<double> perAtomRadii(const core::ParticleData& framework,
                                  const std::vector<double>& perSpeciesRadiusAngstrom) {
    std::vector<double> radii(framework.size());
    for (std::size_t i = 0; i < framework.size(); ++i) {
        const auto species = framework.species[i];
        if (species >= perSpeciesRadiusAngstrom.size()) {
            throw std::invalid_argument(
                "analyzePoreGeometry: framework particle " + std::to_string(i) + " has species " +
                std::to_string(species) + ", but perSpeciesRadiusAngstrom only has " +
                std::to_string(perSpeciesRadiusAngstrom.size()) + " entries");
        }
        radii[i] = perSpeciesRadiusAngstrom[species];
    }
    return radii;
}

/// Converts the Voronoi network's edges into weighted PeriodicEdges. Weight
/// = the edge's own minRadiusAlongEdgeAngstrom (sampled along its actual
/// path at construction time — see VoronoiEdge's doc comment), not just
/// min(endpoint radii): the latter provably misses real constrictions in
/// the interior of anomalously long edges.
std::vector<PeriodicEdge> toWeightedEdges(const VoronoiNetwork& network) {
    std::vector<PeriodicEdge> edges;
    edges.reserve(network.edges.size());
    for (const VoronoiEdge& e : network.edges) {
        edges.push_back(PeriodicEdge{e.from, e.to, e.minRadiusAlongEdgeAngstrom, e.imageOffset});
    }
    return edges;
}

/// Nearest network node (among those with radius >= probeRadiusAngstrom)
/// to `point`, under `lattice`'s minimum-image convention -- the practical
/// point-to-network-region classification behind
/// PoreAnalysis's Monte Carlo accessibility test (see pore_analysis.hpp's
/// doc comment on accessibleSurfaceAreaAngstromSq for why this is a
/// deliberate, documented approximation rather than an exact point-location
/// query against the true Voronoi cell decomposition).
std::size_t nearestQualifyingNode(const VoronoiNetwork& network, const core::Lattice& lattice,
                                   const Vec3& point, double probeRadiusAngstrom) {
    std::size_t best = network.nodes.size();
    double bestDistance = std::numeric_limits<double>::infinity();
    for (std::size_t i = 0; i < network.nodes.size(); ++i) {
        if (network.nodes[i].freeRadiusAngstrom < probeRadiusAngstrom) {
            continue;
        }
        const auto d = lattice.minimumImageDisplacement(point, network.nodes[i].positionAngstrom);
        const double dist = norm(d);
        if (dist < bestDistance) {
            bestDistance = dist;
            best = i;
        }
    }
    return best;
}

} // namespace

PoreGeometry analyzePoreGeometry(const core::ParticleData& framework, const core::Lattice& lattice,
                                  const std::vector<double>& perSpeciesRadiusAngstrom,
                                  const PoreAnalysisOptions& options) {
    if (framework.empty()) {
        throw std::invalid_argument("analyzePoreGeometry: framework must not be empty");
    }
    const std::vector<double> radii = perAtomRadii(framework, perSpeciesRadiusAngstrom);
    const VoronoiNetwork network = buildVoronoiNetwork(framework, lattice, radii);
    const std::vector<PeriodicEdge> weightedEdges = toWeightedEdges(network);

    double maxFreeRadius = network.maxFreeRadiusAlongAnyEdgeAngstrom;
    for (const VoronoiNode& node : network.nodes) {
        maxFreeRadius = std::max(maxFreeRadius, node.freeRadiusAngstrom);
    }
    const double lcd = 2.0 * maxFreeRadius;

    const double pld = 2.0 * percolationThreshold(network.nodes.size(), weightedEdges);

    const std::vector<bool> accessible =
        accessibleComponents(network.nodes.size(), weightedEdges, options.probeRadiusAngstrom);

    // Monte Carlo accessible surface area: sample points on each atom's
    // (radius + probeRadius)-offset sphere, uniformly by direction (Zeo++'s
    // own methodology -- CLAUDE.md section 4's "validated with known
    // deviation" honesty applies here too, see the validation test for the
    // measured tolerance against real Zeo++ output). A sampled point counts
    // toward the area only if a probe centered there doesn't overlap any
    // OTHER atom, and is attributed to accessible or inaccessible surface
    // area by which the nearest qualifying network node belongs to.
    core::Philox4x32Rng rng;
    rng.seed(options.sampleSeed, /*streamIndex=*/0);
    double accessibleArea = 0.0;
    double inaccessibleArea = 0.0;
    for (std::size_t i = 0; i < framework.size(); ++i) {
        const double sphereRadius = radii[i] + options.probeRadiusAngstrom;
        const double atomAreaPerSample =
            4.0 * std::numbers::pi * sphereRadius * sphereRadius /
            static_cast<double>(options.asaSamplesPerAtom);
        const Vec3 center{framework.x[i], framework.y[i], framework.z[i]};
        for (std::size_t s = 0; s < options.asaSamplesPerAtom; ++s) {
            // Uniform point on the unit sphere via the standard
            // z/azimuth parameterization.
            const double z = 2.0 * rng.nextUniform() - 1.0;
            const double phi = 2.0 * std::numbers::pi * rng.nextUniform();
            const double rxy = std::sqrt(std::max(0.0, 1.0 - z * z));
            const Vec3 probeCenter{center[0] + sphereRadius * rxy * std::cos(phi),
                                    center[1] + sphereRadius * rxy * std::sin(phi),
                                    center[2] + sphereRadius * z};

            bool overlaps = false;
            for (std::size_t j = 0; j < framework.size(); ++j) {
                if (j == i) {
                    continue;
                }
                const auto d = lattice.minimumImageDisplacement(
                    probeCenter, {framework.x[j], framework.y[j], framework.z[j]});
                if (norm(d) < radii[j] + options.probeRadiusAngstrom) {
                    overlaps = true;
                    break;
                }
            }
            if (overlaps) {
                continue;
            }

            const std::size_t nearest =
                nearestQualifyingNode(network, lattice, probeCenter, options.probeRadiusAngstrom);
            if (nearest < network.nodes.size() && accessible[nearest]) {
                accessibleArea += atomAreaPerSample;
            } else {
                inaccessibleArea += atomAreaPerSample;
            }
        }
    }

    // Monte Carlo accessible volume: uniform points across the whole unit
    // cell, classified the same way (no atom overlap, then nearest
    // qualifying node's accessibility).
    const double cellVolume = std::abs(lattice.volume());
    const auto& matrix = lattice.matrix();
    const Vec3 a0{matrix[0][0], matrix[0][1], matrix[0][2]};
    const Vec3 a1{matrix[1][0], matrix[1][1], matrix[1][2]};
    const Vec3 a2{matrix[2][0], matrix[2][1], matrix[2][2]};
    std::size_t accessibleCount = 0;
    std::size_t inaccessibleCount = 0;
    for (std::size_t s = 0; s < options.volumeSamplesTotal; ++s) {
        const double u = rng.nextUniform();
        const double v = rng.nextUniform();
        const double w = rng.nextUniform();
        const Vec3 point{u * a0[0] + v * a1[0] + w * a2[0], u * a0[1] + v * a1[1] + w * a2[1],
                          u * a0[2] + v * a1[2] + w * a2[2]};

        bool overlaps = false;
        for (std::size_t j = 0; j < framework.size(); ++j) {
            const auto d =
                lattice.minimumImageDisplacement(point, {framework.x[j], framework.y[j], framework.z[j]});
            if (norm(d) < radii[j] + options.probeRadiusAngstrom) {
                overlaps = true;
                break;
            }
        }
        if (overlaps) {
            continue;
        }

        const std::size_t nearest =
            nearestQualifyingNode(network, lattice, point, options.probeRadiusAngstrom);
        if (nearest < network.nodes.size() && accessible[nearest]) {
            ++accessibleCount;
        } else {
            ++inaccessibleCount;
        }
    }
    const double accessibleVolume =
        cellVolume * static_cast<double>(accessibleCount) / static_cast<double>(options.volumeSamplesTotal);
    const double inaccessibleVolume = cellVolume * static_cast<double>(inaccessibleCount) /
                                       static_cast<double>(options.volumeSamplesTotal);

    return PoreGeometry{lcd,
                         pld,
                         accessibleArea,
                         inaccessibleArea,
                         accessibleVolume,
                         inaccessibleVolume,
                         cellVolume};
}

double zeoPlusPlusDefaultRadiusAngstrom(const std::string& element) {
    // Zeo++ 0.4.7 (lsmo-epfl/zeopp-lsmo) networkinfo.cc,
    // initializeRadTable(), comment "Parameters taken from CCDC tables" --
    // fetched live from the real source (github.com/lsmo-epfl/zeopp-lsmo)
    // during this milestone, not transcribed from memory, so
    // reproductions of Zeo++'s own numbers use the exact same radii.
    static const std::map<std::string, double> table{
        {"H", 1.09},  {"C", 1.70},  {"N", 1.55},  {"O", 1.52},  {"F", 1.47},  {"Al", 2.00},
        {"Si", 2.10}, {"P", 1.80},  {"S", 1.80},  {"Cl", 1.75}, {"Zn", 1.39}, {"Na", 2.27},
        {"Mg", 1.73}, {"K", 2.75},  {"Ca", 2.00},
    };
    const auto it = table.find(element);
    if (it == table.end()) {
        throw std::invalid_argument("zeoPlusPlusDefaultRadiusAngstrom: no entry for element \"" +
                                     element + "\"");
    }
    return it->second;
}

} // namespace aleator::engines
