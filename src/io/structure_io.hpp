#pragma once

#include <filesystem>
#include <string>
#include <vector>

#include "core/geometry/lattice.hpp"
#include "core/math/particle_data.hpp"

namespace aleator::io {

/// Result of reading a structure file: the SoA particle data (CLAUDE.md
/// invariant #2) plus the auxiliary, file-format-specific bookkeeping that
/// doesn't belong in the generic core::ParticleData — which element each
/// species ID refers to, and (for CIF, where this matters) each atom's
/// site occupancy and the asymmetric-unit label it was symmetry-expanded
/// from.
struct StructureData {
    core::ParticleData particles;
    core::Lattice lattice;

    /// speciesSymbols[particles.species[i]] is the element symbol (e.g.
    /// "Si", "O", "Cu") for particle i.
    std::vector<std::string> speciesSymbols;

    /// Parallel to particles: site occupancy (1.0 if the file didn't
    /// specify one). Not resolved into anything physical here — CIF
    /// disorder (multiple atoms at/near the same position, occupancies
    /// summing to ~1) is preserved as-is; deciding what to do with it is a
    /// later milestone's problem.
    std::vector<double> occupancy;

    /// Parallel to particles: the original atom_site label (e.g. "Cu1")
    /// this atom was generated from before symmetry expansion. Distinct
    /// disorder components keep distinct labels, and every symmetry image
    /// of the same asymmetric-unit atom shares one.
    std::vector<std::string> sourceLabel;
};

/// Reads atoms + periodic cell from a CIF file: cell parameters
/// (`_cell_length_*`, `_cell_angle_*`), the symmetry operations
/// (`_symmetry_equiv_pos_as_xyz` or `_space_group_symop_operation_xyz`,
/// defaulting to P1 — identity only — if neither loop is present), and the
/// `_atom_site_*` loop, which is expanded from the asymmetric unit through
/// every symmetry operation with per-site deduplication (so atoms sitting
/// on special positions, e.g. inversion centers, are not double-counted).
///
/// Throws CifParseError on any structurally malformed input, and
/// std::filesystem::filesystem_error if `file` cannot be opened. Never a
/// crash, and never a silently wrong parse (CLAUDE.md invariant #7).
[[nodiscard]] StructureData readCif(const std::filesystem::path& file);

/// Writes atoms + periodic cell to a PDB file (CRYST1 record for the cell).
/// Declared only — writing is out of scope for this milestone. Throws
/// NotImplemented.
void writePdb(const std::filesystem::path& file, const core::ParticleData& particles,
              const core::Lattice& lattice);

/// Writes atoms + periodic cell in LAMMPS `data` file format. Declared
/// only. Throws NotImplemented.
void writeLammpsData(const std::filesystem::path& file, const core::ParticleData& particles,
                      const core::Lattice& lattice);

/// Appends one frame (positions only) to a trajectory file. Declared only.
/// Throws NotImplemented.
void appendTrajectoryFrame(const std::filesystem::path& file, const core::ParticleData& particles,
                            const core::Lattice& lattice);

} // namespace aleator::io
