#include "core/neighbor/verlet_list.hpp"

#include <array>
#include <stdexcept>

#include "core/neighbor/cell_list.hpp"

namespace aleator::core {

namespace {

double squaredDistance(const Lattice& lattice, const std::array<double, 3>& a,
                        const std::array<double, 3>& b) {
    const auto d = lattice.minimumImageDisplacement(a, b);
    return d[0] * d[0] + d[1] * d[1] + d[2] * d[2];
}

} // namespace

VerletList::VerletList(double skinDistance) : skin_(skinDistance) {
    if (!(skinDistance > 0.0)) {
        throw std::invalid_argument("VerletList: skinDistance must be positive");
    }
}

void VerletList::build(const ParticleData& particles, const Lattice& lattice, double cutoff) {
    if (!(cutoff > 0.0)) {
        throw std::invalid_argument("VerletList::build: cutoff must be positive");
    }
    // Rejects cutoff > L_perp/2 with a clear, direction-naming error
    // (CLAUDE.md #10) before anything else happens.
    lattice.validateCutoff(cutoff);

    cutoff_ = cutoff;
    lattice_ = lattice;
    particles_ = &particles;

    referenceX_.assign(particles.x.begin(), particles.x.end());
    referenceY_.assign(particles.y.begin(), particles.y.end());
    referenceZ_.assign(particles.z.begin(), particles.z.end());

    CellList cellList;
    cellList.build(particles, lattice, cutoff + skin_);

    const double candidateCutoffSquared = (cutoff + skin_) * (cutoff + skin_);
    candidatePairs_.clear();
    cellList.forEachCandidatePair([&](std::size_t i, std::size_t j) {
        const std::array<double, 3> pi{particles.x[i], particles.y[i], particles.z[i]};
        const std::array<double, 3> pj{particles.x[j], particles.y[j], particles.z[j]};
        if (squaredDistance(lattice, pi, pj) <= candidateCutoffSquared) {
            candidatePairs_.push_back({i, j});
        }
    });

    built_ = true;
}

bool VerletList::needsRebuild(const ParticleData& particles) const {
    if (!built_) {
        return true;
    }
    if (particles.size() != referenceX_.size()) {
        // Particle count changed (e.g. GCMC insertion/deletion) since the
        // last build(): the candidate superset no longer even indexes the
        // right set of particles, so a rebuild is unconditionally required.
        return true;
    }

    const double halfSkinSquared = (skin_ * 0.5) * (skin_ * 0.5);
    for (std::size_t i = 0; i < particles.size(); ++i) {
        const std::array<double, 3> reference{referenceX_[i], referenceY_[i], referenceZ_[i]};
        const std::array<double, 3> current{particles.x[i], particles.y[i], particles.z[i]};
        if (squaredDistance(lattice_, reference, current) >= halfSkinSquared) {
            return true;
        }
    }
    return false;
}

const std::vector<NeighborPair>& VerletList::pairs() const {
    exactPairsCache_.clear();
    if (particles_ == nullptr) {
        return exactPairsCache_;
    }
    const double cutoffSquared = cutoff_ * cutoff_;
    for (const auto& candidate : candidatePairs_) {
        const std::array<double, 3> pi{particles_->x[candidate.i], particles_->y[candidate.i],
                                        particles_->z[candidate.i]};
        const std::array<double, 3> pj{particles_->x[candidate.j], particles_->y[candidate.j],
                                        particles_->z[candidate.j]};
        if (squaredDistance(lattice_, pi, pj) <= cutoffSquared) {
            exactPairsCache_.push_back(candidate);
        }
    }
    return exactPairsCache_;
}

} // namespace aleator::core
