#pragma once

#include <array>

#include "core/geometry/lattice.hpp"
#include "core/math/particle_data.hpp"

namespace aleator::core {

/// Minimum (nx, ny, nz), each >= 1, such that replicating `lattice` by this
/// amount along each axis makes every perpendicular width strictly exceed
/// `2 * cutoffAngstrom` -- the minimum-image safety margin
/// `Lattice::validateCutoff` already enforces elsewhere, computed per axis
/// from `Lattice::perpendicularWidth` (NOT raw lattice vector length: for a
/// triclinic cell these differ, and using `|a_i|` is exactly the kind of
/// "assumes orthorhombic" bug CLAUDE.md invariant 6 warns about). Each axis
/// is independent: replicating only along axis `i` scales `perpendicularWidth(i)`
/// by exactly that axis's replication factor and leaves the other two
/// unchanged (the cell volume scales linearly in one lattice vector while
/// the cross product defining the other two axes' perpendicular width is
/// untouched), so there is no coupling between axes to solve for jointly.
/// Throws std::invalid_argument if `cutoffAngstrom` is not positive.
[[nodiscard]] std::array<int, 3> minimumSupercellReplication(const Lattice& lattice,
                                                               double cutoffAngstrom);

/// A replicated system: `particles` under `lattice`.
struct SupercellResult {
    ParticleData particles;
    Lattice lattice;
};

/// Replicates `particles` under `lattice` by (nx, ny, nz), each >= 1 --
/// an exact operation (integer lattice-vector translations of the original
/// Cartesian positions, not a fractional round-trip: image (i,j,k)'s copy
/// of particle p sits at p's original position plus i*a0 + j*a1 + k*a2,
/// using `lattice`'s own row vectors), returning the replicated particles
/// and the correspondingly-scaled lattice (row i scaled by n_i). (1,1,1)
/// is a true no-op: every translation offset is exactly 0.0, so particle
/// positions and the lattice matrix come out bit-identical to the input,
/// not merely close -- this is deliberate, since CLAUDE.md section 5's
/// exact-refactor bar for anything meant to be exact requires it, and this
/// specific identity case is what lets existing single-cell examples/tests
/// keep reproducing their exact numbers when the CLI's automatic
/// replication computes (1,1,1) for them.
///
/// Throws std::invalid_argument if any of nx, ny, nz is not >= 1.
[[nodiscard]] SupercellResult replicateSupercell(const ParticleData& particles,
                                                   const Lattice& lattice, int nx, int ny, int nz);

} // namespace aleator::core
