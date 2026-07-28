#pragma once

#include <cstdint>
#include <vector>

namespace aleator::engines {

/// One interaction site of a rigid adsorbate molecule template, in the
/// molecule's own local (unrotated) frame — offset from the molecule's
/// reference point (Å), the LJ species index it maps to in the
/// ForceField's per-species parameter table, its point charge (e), and
/// mass (amu). A single-site species (e.g. united-atom CH4) has exactly
/// one site at the origin.
struct MoleculeSite {
    double x = 0.0;
    double y = 0.0;
    double z = 0.0;
    std::uint32_t ljSpecies = 0;
    double charge = 0.0;
    double mass = 0.0;
};

/// A rigid multi-site adsorbate template. GCMC insertion places a copy at
/// a random position and (for more than one site) random orientation;
/// translation/rotation moves rigidly move all of an existing instance's
/// sites together.
struct MoleculeSpecies {
    std::vector<MoleculeSite> sites;
};

} // namespace aleator::engines
