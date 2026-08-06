#include "engines/monte_carlo/framework_energy_grid.hpp"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstring>
#include <fstream>
#include <sstream>
#include <stdexcept>
#include <type_traits>

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

// ---- raw binary I/O helpers -------------------------------------------
// Deliberately plain, fixed-width, no external serialization dependency:
// the cache format only ever needs to round-trip within this codebase's
// own reader/writer, both compiled from the same source, on the portability
// targets CLAUDE.md section 6 lists (all little-endian in practice; a
// byte-order mismatch is out of scope the same way it already is for every
// other raw-binary artifact in this codebase).
template <typename T>
void writePod(std::ostream& out, const T& value) {
    static_assert(std::is_trivially_copyable_v<T>);
    out.write(reinterpret_cast<const char*>(&value), sizeof(T));
}
template <typename T>
void readPod(std::istream& in, T& value) {
    static_assert(std::is_trivially_copyable_v<T>);
    in.read(reinterpret_cast<char*>(&value), sizeof(T));
}

void writeDoubleVector(std::ostream& out, const std::vector<double>& v) {
    writePod(out, static_cast<std::uint64_t>(v.size()));
    if (!v.empty()) {
        out.write(reinterpret_cast<const char*>(v.data()),
                   static_cast<std::streamsize>(v.size() * sizeof(double)));
    }
}
std::vector<double> readDoubleVector(std::istream& in) {
    std::uint64_t n = 0;
    readPod(in, n);
    std::vector<double> v(n);
    if (n > 0) {
        in.read(reinterpret_cast<char*>(v.data()), static_cast<std::streamsize>(n * sizeof(double)));
    }
    return v;
}

void writeU32Vector(std::ostream& out, const std::vector<std::uint32_t>& v) {
    writePod(out, static_cast<std::uint64_t>(v.size()));
    if (!v.empty()) {
        out.write(reinterpret_cast<const char*>(v.data()),
                   static_cast<std::streamsize>(v.size() * sizeof(std::uint32_t)));
    }
}
std::vector<std::uint32_t> readU32Vector(std::istream& in) {
    std::uint64_t n = 0;
    readPod(in, n);
    std::vector<std::uint32_t> v(n);
    if (n > 0) {
        in.read(reinterpret_cast<char*>(v.data()),
                static_cast<std::streamsize>(n * sizeof(std::uint32_t)));
    }
    return v;
}

void writeKey(std::ostream& out, const FrameworkEnergyGridCacheKey& key) {
    writeDoubleVector(out, key.atomX);
    writeDoubleVector(out, key.atomY);
    writeDoubleVector(out, key.atomZ);
    writeU32Vector(out, key.atomSpecies);
    for (const auto& row : key.latticeMatrix) {
        for (double v : row) {
            writePod(out, v);
        }
    }
    writePod(out, static_cast<std::uint64_t>(key.perSpeciesLJ.size()));
    for (const auto& p : key.perSpeciesLJ) {
        writePod(out, p.epsilon);
        writePod(out, p.sigma);
    }
    writePod(out, key.cutoffAngstrom);
    writePod(out, key.truncation);
    writeU32Vector(out, key.guestSpeciesIds);
    writePod(out, key.spacingAngstrom);
    writePod(out, key.energyCapKelvin);
}

FrameworkEnergyGridCacheKey readKey(std::istream& in) {
    FrameworkEnergyGridCacheKey key;
    key.atomX = readDoubleVector(in);
    key.atomY = readDoubleVector(in);
    key.atomZ = readDoubleVector(in);
    key.atomSpecies = readU32Vector(in);
    for (auto& row : key.latticeMatrix) {
        for (double& v : row) {
            readPod(in, v);
        }
    }
    std::uint64_t perSpeciesCount = 0;
    readPod(in, perSpeciesCount);
    key.perSpeciesLJ.resize(perSpeciesCount);
    for (auto& p : key.perSpeciesLJ) {
        readPod(in, p.epsilon);
        readPod(in, p.sigma);
    }
    readPod(in, key.cutoffAngstrom);
    readPod(in, key.truncation);
    key.guestSpeciesIds = readU32Vector(in);
    readPod(in, key.spacingAngstrom);
    readPod(in, key.energyCapKelvin);
    return key;
}

// FNV-1a 64-bit -- see FrameworkEnergyGridCacheKey::hash()'s doc comment
// for why this is only ever used to pick a filename, never as the sole
// correctness check.
std::uint64_t fnv1a64(const void* data, std::size_t size) {
    constexpr std::uint64_t kOffsetBasis = 0xcbf29ce484222325ULL;
    constexpr std::uint64_t kPrime = 0x100000001b3ULL;
    std::uint64_t hash = kOffsetBasis;
    const auto* bytes = static_cast<const unsigned char*>(data);
    for (std::size_t i = 0; i < size; ++i) {
        hash ^= static_cast<std::uint64_t>(bytes[i]);
        hash *= kPrime;
    }
    return hash;
}

/// Reads and returns (formatVersion, key) from `path`'s header, without
/// touching the grid table data that follows. Throws std::runtime_error if
/// the file can't be opened or is truncated mid-header (a corrupt file is
/// exactly the kind of thing that must hard-error, not be treated as "no
/// cache entry").
std::pair<std::uint32_t, FrameworkEnergyGridCacheKey> readHeader(const std::filesystem::path& path) {
    std::ifstream in(path, std::ios::binary);
    if (!in) {
        throw std::runtime_error("FrameworkEnergyGrid cache: could not open " + path.string());
    }
    std::uint32_t version = 0;
    readPod(in, version);
    FrameworkEnergyGridCacheKey key = readKey(in);
    if (!in) {
        throw std::runtime_error("FrameworkEnergyGrid cache: " + path.string() +
                                  " is truncated or corrupt (header read failed)");
    }
    return {version, std::move(key)};
}

} // namespace

std::uint64_t FrameworkEnergyGridCacheKey::hash() const {
    std::ostringstream oss(std::ios::binary);
    writeKey(oss, *this);
    const std::string bytes = oss.str();
    return fnv1a64(bytes.data(), bytes.size());
}

FrameworkEnergyGridCacheKey makeFrameworkEnergyGridCacheKey(
    const forcefield::LennardJones& forceField, const core::ParticleData& frameworkParticles,
    const core::Lattice& lattice, const std::vector<std::uint32_t>& guestSpeciesIds,
    double spacingAngstrom) {
    FrameworkEnergyGridCacheKey key;
    const std::size_t n = frameworkParticles.size();
    key.atomX.assign(frameworkParticles.x.begin(), frameworkParticles.x.end());
    key.atomY.assign(frameworkParticles.y.begin(), frameworkParticles.y.end());
    key.atomZ.assign(frameworkParticles.z.begin(), frameworkParticles.z.end());
    key.atomSpecies.assign(frameworkParticles.species.begin(), frameworkParticles.species.end());
    (void)n;
    key.latticeMatrix = lattice.matrix();
    key.perSpeciesLJ = forceField.perSpecies();
    key.cutoffAngstrom = forceField.cutoff();
    key.truncation = static_cast<std::uint8_t>(forceField.truncation());
    key.guestSpeciesIds = guestSpeciesIds;
    std::sort(key.guestSpeciesIds.begin(), key.guestSpeciesIds.end());
    key.guestSpeciesIds.erase(std::unique(key.guestSpeciesIds.begin(), key.guestSpeciesIds.end()),
                               key.guestSpeciesIds.end());
    key.spacingAngstrom = spacingAngstrom;
    key.energyCapKelvin = FrameworkEnergyGrid::kEnergyCapKelvin;
    return key;
}

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

void FrameworkEnergyGrid::serialize(const std::filesystem::path& path,
                                     const FrameworkEnergyGridCacheKey& key) const {
    std::filesystem::create_directories(path.parent_path());
    std::ofstream out(path, std::ios::binary | std::ios::trunc);
    if (!out) {
        throw std::runtime_error("FrameworkEnergyGrid::serialize: could not open " + path.string() +
                                  " for writing");
    }
    writePod(out, kFrameworkEnergyGridCacheFormatVersion);
    writeKey(out, key);
    for (int axis = 0; axis < 3; ++axis) {
        writePod(out, pointCounts_[static_cast<std::size_t>(axis)]);
    }
    writePod(out, static_cast<std::uint64_t>(tables_.size()));
    // key.guestSpeciesIds is the sorted/deduplicated set makeFrameworkEnergyGridCacheKey
    // produced; tables_ was built from that same set (see the constructor),
    // so iterating tables_ directly (rather than re-deriving from
    // key.guestSpeciesIds) keeps this the single source of truth for "which
    // tables actually exist" and can't drift from what was really built.
    for (const auto& [speciesId, table] : tables_) {
        writePod(out, speciesId);
        writeDoubleVector(out, table.values);
    }
    if (!out) {
        throw std::runtime_error("FrameworkEnergyGrid::serialize: write failed for " + path.string());
    }
}

std::filesystem::path FrameworkEnergyGrid::cacheFilePath(const std::filesystem::path& cacheDirectory,
                                                           const FrameworkEnergyGridCacheKey& key) {
    char hex[17];
    std::snprintf(hex, sizeof(hex), "%016llx", static_cast<unsigned long long>(key.hash()));
    return cacheDirectory /
           ("framework_energy_grid_v" + std::to_string(kFrameworkEnergyGridCacheFormatVersion) + "_" +
            std::string(hex) + ".bin");
}

bool FrameworkEnergyGrid::cacheHit(const std::filesystem::path& cacheDirectory,
                                    const FrameworkEnergyGridCacheKey& key) {
    const auto path = cacheFilePath(cacheDirectory, key);
    if (!std::filesystem::exists(path)) {
        return false;
    }
    const auto [version, storedKey] = readHeader(path);
    if (version != kFrameworkEnergyGridCacheFormatVersion) {
        throw std::runtime_error(
            "FrameworkEnergyGrid cache: " + path.string() + " has format version " +
            std::to_string(version) + ", this build expects " +
            std::to_string(kFrameworkEnergyGridCacheFormatVersion) +
            " -- refusing to reinterpret it. Delete the file (or the cache directory) to force a "
            "rebuild.");
    }
    if (!(storedKey == key)) {
        throw std::runtime_error(
            "FrameworkEnergyGrid cache: " + path.string() +
            " exists and has a matching content-hash filename, but its stored key does not "
            "exactly match the current inputs (framework/lattice/force-field/spacing) -- this "
            "is either a hash collision or a corrupted/foreign cache file, and refusing to "
            "silently reuse or silently rebuild over it is the only safe response. Delete the "
            "file if you are sure it is stale.");
    }
    return true;
}

std::optional<FrameworkEnergyGrid> FrameworkEnergyGrid::loadFromCache(
    const std::filesystem::path& cacheDirectory, const FrameworkEnergyGridCacheKey& key) {
    const auto path = cacheFilePath(cacheDirectory, key);
    if (!std::filesystem::exists(path)) {
        return std::nullopt;
    }

    std::ifstream in(path, std::ios::binary);
    if (!in) {
        throw std::runtime_error("FrameworkEnergyGrid::loadFromCache: could not open " +
                                  path.string());
    }
    std::uint32_t version = 0;
    readPod(in, version);
    if (version != kFrameworkEnergyGridCacheFormatVersion) {
        throw std::runtime_error(
            "FrameworkEnergyGrid cache: " + path.string() + " has format version " +
            std::to_string(version) + ", this build expects " +
            std::to_string(kFrameworkEnergyGridCacheFormatVersion) +
            " -- refusing to reinterpret it. Delete the file (or the cache directory) to force a "
            "rebuild.");
    }
    FrameworkEnergyGridCacheKey storedKey = readKey(in);
    if (!in) {
        throw std::runtime_error("FrameworkEnergyGrid cache: " + path.string() +
                                  " is truncated or corrupt (key read failed)");
    }
    if (!(storedKey == key)) {
        throw std::runtime_error(
            "FrameworkEnergyGrid cache: " + path.string() +
            " exists and has a matching content-hash filename, but its stored key does not "
            "exactly match the current inputs -- refusing to silently reuse a possibly-stale or "
            "colliding cache entry. Delete the file if you are sure it is stale.");
    }

    FrameworkEnergyGrid grid;
    grid.lattice_ = core::Lattice(key.latticeMatrix);
    grid.spacingAngstrom_ = key.spacingAngstrom;
    for (int axis = 0; axis < 3; ++axis) {
        readPod(in, grid.pointCounts_[static_cast<std::size_t>(axis)]);
    }
    std::uint64_t tableCount = 0;
    readPod(in, tableCount);
    for (std::uint64_t i = 0; i < tableCount; ++i) {
        std::uint32_t speciesId = 0;
        readPod(in, speciesId);
        Table table;
        table.values = readDoubleVector(in);
        grid.tables_.emplace(speciesId, std::move(table));
    }
    if (!in) {
        throw std::runtime_error("FrameworkEnergyGrid cache: " + path.string() +
                                  " is truncated or corrupt (grid data read failed)");
    }
    return grid;
}

FrameworkEnergyGrid loadOrBuildFrameworkEnergyGrid(const forcefield::LennardJones& forceField,
                                                    const core::ParticleData& frameworkParticles,
                                                    const core::Lattice& lattice,
                                                    const FrameworkEnergyGridCacheRequest& request,
                                                    bool* wasHit) {
    const FrameworkEnergyGridCacheKey key = makeFrameworkEnergyGridCacheKey(
        forceField, frameworkParticles, lattice, request.guestSpeciesIds, request.spacingAngstrom);

    auto cached = FrameworkEnergyGrid::loadFromCache(request.cacheDirectory, key);
    if (cached.has_value()) {
        if (wasHit != nullptr) {
            *wasHit = true;
        }
        return std::move(*cached);
    }

    if (wasHit != nullptr) {
        *wasHit = false;
    }
    FrameworkEnergyGrid grid(forceField, frameworkParticles, lattice, request.guestSpeciesIds,
                              request.spacingAngstrom);
    const auto path = FrameworkEnergyGrid::cacheFilePath(request.cacheDirectory, key);
    grid.serialize(path, key);
    return grid;
}

} // namespace aleator::engines
