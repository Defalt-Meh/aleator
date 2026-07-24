#pragma once

#include <cstdint>
#include <filesystem>
#include <string>
#include <unordered_map>

#include "forcefield/pairwise/lennard_jones.hpp"

namespace aleator::forcefield {

/// Force-field parameter set as loaded from a parameter file: per-species
/// Lennard-Jones parameters keyed by species name, plus whatever additional
/// sections (electrostatics, bonded terms) get added as those force fields
/// are implemented.
struct ParameterSet {
    std::unordered_map<std::string, LennardJonesParameters> lennardJones;
};

/// Loads a force-field parameter file (TOML; see io/ for the general TOML
/// config reader). Declared only: parsing is straightforward, but silently
/// mis-mapping a column (e.g. swapping epsilon/sigma units) is exactly the
/// kind of "smooth, plausible, wrong" failure CLAUDE.md warns about, so this
/// waits for a round-trip test against a known parameter file before it is
/// implemented. Throws NotImplemented.
[[nodiscard]] ParameterSet loadParameterSet(const std::filesystem::path& file);

} // namespace aleator::forcefield
