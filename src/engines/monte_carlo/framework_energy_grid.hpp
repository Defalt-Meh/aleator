#pragma once

#include <array>
#include <cstddef>
#include <cstdint>
#include <filesystem>
#include <optional>
#include <unordered_map>
#include <vector>

#include "core/geometry/lattice.hpp"
#include "core/math/particle_data.hpp"
#include "forcefield/pairwise/lennard_jones.hpp"

namespace aleator::engines {

/// On-disk cache format version. Bumped whenever the serialized layout, the
/// key hashing/comparison scheme, or the grid-construction algorithm itself
/// changes in a way that would make an old cache file's bytes mean
/// something different. A version mismatch is a hard error on load (see
/// FrameworkEnergyGrid::loadFromCache) -- never silently reinterpreted,
/// never silently rebuilt as if nothing were wrong.
inline constexpr std::uint32_t kFrameworkEnergyGridCacheFormatVersion = 1;

/// Every input that determines a FrameworkEnergyGrid's contents, bundled so
/// it can be (a) hashed, for a content-addressed cache filename, and (b)
/// stored verbatim in the cache file and compared field-by-field on load.
/// The hash alone is never trusted for correctness -- only for choosing a
/// filename -- because a hash collision (or a future hashing bug) silently
/// reusing the wrong file's data would be exactly the "stale grid silently
/// producing wrong energies" failure CLAUDE.md section 5 calls the worst
/// failure mode in this project. The exact field-by-field `operator==`
/// below is what actually gates reuse.
///
/// Every field here is a plain, exactly-comparable value (no floating-point
/// arithmetic happens between "compute the key" and "compare two keys" --
/// they are the same literal inputs a caller already has, e.g. the CIF's
/// own parsed atom positions), so comparing doubles with `==` here is
/// correct and intentional, not a numerical-tolerance question.
struct FrameworkEnergyGridCacheKey {
    std::vector<double> atomX, atomY, atomZ;
    std::vector<std::uint32_t> atomSpecies;
    std::array<std::array<double, 3>, 3> latticeMatrix{};
    std::vector<forcefield::LennardJonesParameters> perSpeciesLJ;
    double cutoffAngstrom = 0.0;
    std::uint8_t truncation = 0;
    std::vector<std::uint32_t> guestSpeciesIds; // sorted + deduplicated
    double spacingAngstrom = 0.0;
    double energyCapKelvin = 0.0;

    [[nodiscard]] bool operator==(const FrameworkEnergyGridCacheKey&) const = default;

    /// FNV-1a 64-bit over a canonical byte serialization of every field
    /// above, in declaration order. Deterministic across runs/platforms
    /// with the same bit-representation of double (IEEE754, universal on
    /// every CLAUDE.md section 6 portability target) -- used ONLY to pick
    /// a cache filename (content-addressed storage), never as the sole
    /// correctness check (see the struct's own doc comment above).
    [[nodiscard]] std::uint64_t hash() const;
};

[[nodiscard]] FrameworkEnergyGridCacheKey makeFrameworkEnergyGridCacheKey(
    const forcefield::LennardJones& forceField, const core::ParticleData& frameworkParticles,
    const core::Lattice& lattice, const std::vector<std::uint32_t>& guestSpeciesIds,
    double spacingAngstrom);

/// Precomputed guest-host Lennard-Jones energy field for a rigid framework
/// (CLAUDE.md section 5 performance milestone): "the guest-host energy is a
/// fixed scalar field that can be tabulated once and interpolated... the
/// single largest available speedup" for GCMC against a host that never
/// moves. One table per guest LJ species (a multi-site adsorbate like CO2
/// needs one per site type), each a regular periodic grid of energies —
/// sum over every framework atom of pairEnergy(distance, mixedParameters)
/// — evaluated at that grid node, trilinearly interpolated at arbitrary
/// query positions in between.
///
/// Construction cost is O(gridPointCount * frameworkAtomCount) per guest
/// species — real work, meant to be paid once and reused: in-process
/// across several MonteCarloEngines built against the same
/// framework/force field (e.g. an isotherm's pressure points), and now
/// ACROSS PROCESSES too, via serialize()/loadFromCache() below -- a grid's
/// entire purpose is amortization, which a build that can never persist
/// can never actually achieve regardless of how many pressure points
/// follow in a single run (see CLAUDE.md section 0/5 for the earlier
/// session's measurement that first surfaced this gap).
///
/// Not an exact refactor (unlike core::CellList's guest-guest wiring):
/// interpolated energies differ from a direct evaluation at the same point
/// by a real, measured, documented amount that shrinks as spacingAngstrom()
/// shrinks — see the accuracy-vs-spacing measurement in
/// tests/validation/test_framework_energy_grid.cc and the default chosen in
/// engines/monte_carlo/framework_energy_grid.cc's kDefaultSpacingAngstrom
/// doc comment.
class FrameworkEnergyGrid {
public:
    /// Every stored grid node is clamped to at most this (K) before being
    /// interpolated -- see framework_energy_grid.cc's doc comment: raw LJ
    /// repulsion at a node very close to a framework atom's core can
    /// legitimately compute 1e17+ K, and interpolating that against normal
    /// neighboring nodes produces physically meaningless multi-order-of-
    /// magnitude output. At simulation temperatures (~300 K),
    /// exp(-kEnergyCapKelvin / T) is already indistinguishable from exactly
    /// zero acceptance probability, so this cannot change which trial moves
    /// get accepted -- only avoids storing/interpolating meaningless raw
    /// magnitudes. Public so callers comparing interpolated energies
    /// against a direct evaluation (e.g. this class's own accuracy test)
    /// can apply the same cap for a physically-meaningful comparison,
    /// rather than being dominated by cosmetically enormous but Boltzmann-
    /// irrelevant deviations deep in the repulsive core.
    static constexpr double kEnergyCapKelvin = 1.0e5;

    /// Builds one grid per species id in `guestSpeciesIds` (duplicates
    /// ignored). `frameworkParticles` must contain only the fixed host
    /// atoms (never adsorbate sites) -- this grid has no notion of guest
    /// occupancy, by design (a rigid framework's guest-host field doesn't
    /// depend on which/how-many guests are currently present). Throws
    /// std::invalid_argument if spacingAngstrom is not positive.
    FrameworkEnergyGrid(const forcefield::LennardJones& forceField,
                         const core::ParticleData& frameworkParticles, const core::Lattice& lattice,
                         const std::vector<std::uint32_t>& guestSpeciesIds, double spacingAngstrom);

    /// Trilinearly-interpolated guest-host energy (K) for a guest of
    /// `guestSpeciesId` at `position` (Cartesian Å, any periodic image).
    /// O(1): 8 grid-node lookups plus a weighted sum, replacing the
    /// O(frameworkAtomCount) direct scan this class exists to amortize.
    /// Throws std::out_of_range if `guestSpeciesId` wasn't included in
    /// `guestSpeciesIds` at construction.
    [[nodiscard]] double interpolate(std::uint32_t guestSpeciesId,
                                      const std::array<double, 3>& position) const;

    [[nodiscard]] double spacingAngstrom() const noexcept { return spacingAngstrom_; }
    [[nodiscard]] std::array<int, 3> pointCounts() const noexcept { return pointCounts_; }

    /// Writes this grid, plus `key` (the cache key that determines it --
    /// see FrameworkEnergyGridCacheKey; not reconstructed from this object,
    /// since a FrameworkEnergyGrid doesn't retain the framework/force-field
    /// inputs it was built from), to `path` -- creates parent directories
    /// if needed. Overwrites `path` if it already exists. Callers normally
    /// don't call this directly; loadOrBuildFrameworkEnergyGrid() does, on
    /// a cache miss, with the same key it used for the miss check.
    void serialize(const std::filesystem::path& path, const FrameworkEnergyGridCacheKey& key) const;

    /// Deterministic, content-addressed filename for `key` inside
    /// `cacheDirectory` -- does not touch the filesystem. Two different
    /// keys (differing in even one field) are, for all practical purposes,
    /// guaranteed a different path (64-bit FNV-1a hash); a version bump to
    /// kFrameworkEnergyGridCacheFormatVersion is folded into the name too,
    /// so an old-format file can never collide with a new-format one.
    [[nodiscard]] static std::filesystem::path cacheFilePath(
        const std::filesystem::path& cacheDirectory, const FrameworkEnergyGridCacheKey& key);

    /// Cheap hit/miss check for `key` inside `cacheDirectory`, reading only
    /// the file's header (not the, possibly large, grid table data) --
    /// meant for `--dry-run`, which must report whether a real run would
    /// hit or build without paying either the build cost or the full
    /// deserialization cost. Returns false for a normal cache miss (no
    /// file at cacheFilePath(cacheDirectory, key)). Throws
    /// std::runtime_error if a file exists there but its stored format
    /// version or key does not exactly match `key` -- the same "hard
    /// error, never silently reinterpreted" rule loadFromCache() applies,
    /// since silently reporting a cache miss here (when a file exists but
    /// doesn't verify) could lead a caller to instead build fresh AND then
    /// silently overwrite what might be a legitimately-different cached
    /// grid some other job still depends on.
    [[nodiscard]] static bool cacheHit(const std::filesystem::path& cacheDirectory,
                                        const FrameworkEnergyGridCacheKey& key);

    /// Loads a previously-serialized grid from `cacheFilePath(cacheDirectory,
    /// key)`, verifying the on-disk format version and the stored cache key
    /// both exactly match `key` before returning it. Returns std::nullopt
    /// for a normal cache miss (no file at that path). Throws
    /// std::runtime_error -- a hard error, never a silent rebuild and never
    /// a silent reuse -- if a file exists there but its version or key does
    /// not match: this indicates a hash collision against a different key,
    /// a corrupted file, or a tampered/foreign cache directory, none of
    /// which may ever be silently papered over (CLAUDE.md section 5: "a
    /// stale grid silently producing wrong energies is the worst failure
    /// mode in this project").
    [[nodiscard]] static std::optional<FrameworkEnergyGrid> loadFromCache(
        const std::filesystem::path& cacheDirectory, const FrameworkEnergyGridCacheKey& key);

private:
    FrameworkEnergyGrid() = default; // loadFromCache() only

    struct Table {
        std::vector<double> values; // flattened, size nx*ny*nz, index = (iz*ny + iy)*nx + ix
    };

    [[nodiscard]] std::size_t flatten(int ix, int iy, int iz) const noexcept {
        return (static_cast<std::size_t>(iz) * static_cast<std::size_t>(pointCounts_[1]) +
                static_cast<std::size_t>(iy)) *
                   static_cast<std::size_t>(pointCounts_[0]) +
               static_cast<std::size_t>(ix);
    }

    core::Lattice lattice_;
    double spacingAngstrom_ = 0.0;
    std::array<int, 3> pointCounts_{1, 1, 1};
    std::unordered_map<std::uint32_t, Table> tables_;
};

/// Cache directory + guest species ids to build/load for, bundled since
/// every call site (CLI dry-run reporting, CLI real-run build, benchmarks)
/// needs the same tuple.
struct FrameworkEnergyGridCacheRequest {
    std::filesystem::path cacheDirectory;
    std::vector<std::uint32_t> guestSpeciesIds;
    double spacingAngstrom = 0.0;
};

/// Orchestrates the cache: computes the key, checks
/// FrameworkEnergyGrid::loadFromCache() first, and only builds from scratch
/// (then serializes the result for next time) on a genuine miss. `wasHit`,
/// if non-null, is set to whether the cache was hit -- for real-run
/// logging and for tests distinguishing the two paths without re-deriving
/// the key themselves.
[[nodiscard]] FrameworkEnergyGrid loadOrBuildFrameworkEnergyGrid(
    const forcefield::LennardJones& forceField, const core::ParticleData& frameworkParticles,
    const core::Lattice& lattice, const FrameworkEnergyGridCacheRequest& request,
    bool* wasHit = nullptr);

} // namespace aleator::engines
