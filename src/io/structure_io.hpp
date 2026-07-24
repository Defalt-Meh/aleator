#pragma once

#include <filesystem>
#include <utility>

#include "core/geometry/lattice.hpp"
#include "core/math/particle_data.hpp"

namespace aleator::io {

/// Framework/structure readers and writers. Declared only: file-format
/// parsing is a correctness boundary (CLAUDE.md #2.4 — "every I/O boundary
/// converts explicitly"; a wrong CIF symmetry expansion or a silently
/// mis-parsed PDB column is exactly the class of bug the "no stubs
/// presented as complete" invariant is meant to prevent). Each throws
/// NotImplemented.

/// Reads atoms + periodic cell from a CIF file.
[[nodiscard]] std::pair<core::ParticleData, core::Lattice> readCif(
    const std::filesystem::path& file);

/// Writes atoms + periodic cell to a PDB file (CRYST1 record for the cell).
void writePdb(const std::filesystem::path& file, const core::ParticleData& particles,
              const core::Lattice& lattice);

/// Writes atoms + periodic cell in LAMMPS `data` file format.
void writeLammpsData(const std::filesystem::path& file, const core::ParticleData& particles,
                      const core::Lattice& lattice);

/// Appends one frame (positions only) to a trajectory file.
void appendTrajectoryFrame(const std::filesystem::path& file, const core::ParticleData& particles,
                            const core::Lattice& lattice);

} // namespace aleator::io
