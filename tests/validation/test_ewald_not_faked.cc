// Physics validation tests live under this label and are gated on
// literature/analytic reference values (CLAUDE.md #4) — none exist yet
// because no physics is implemented in this build.
//
// This test still belongs here rather than under [unit]: it guards
// invariant #7 ("no stubs presented as complete... never a hardcoded
// plausible return value") specifically for Ewald, the component CLAUDE.md
// #4 calls out as highest-risk precisely because a subtly wrong
// implementation produces smooth, plausible, wrong isotherms. Until the
// Madelung-constant validation test (NaCl rocksalt, 1.747564594633...,
// 1e-6 relative) is wired up and passing, computeEnergy must throw rather
// than return anything that could be mistaken for a real result.

#include <catch2/catch_test_macros.hpp>

#include <vector>

#include "core/exceptions.hpp"
#include "core/geometry/lattice.hpp"
#include "core/math/particle_data.hpp"
#include "forcefield/electrostatics/ewald.hpp"

namespace {

/// Minimal NeighborList that reports zero pairs — sufficient to reach
/// Ewald::computeEnergy without needing a real neighbor list yet.
class EmptyNeighborList final : public aleator::core::NeighborList {
public:
    void build(const aleator::core::ParticleData&, const aleator::core::Lattice&,
               double) override {}
    [[nodiscard]] bool needsRebuild(const aleator::core::ParticleData&) const override {
        return false;
    }
    [[nodiscard]] const std::vector<aleator::core::NeighborPair>& pairs() const override {
        return pairs_;
    }

private:
    std::vector<aleator::core::NeighborPair> pairs_;
};

} // namespace

TEST_CASE("Ewald::computeEnergy throws rather than returning a fabricated energy",
          "[validation][electrostatics]") {
    const aleator::forcefield::Ewald ewald(/*alpha=*/0.3, /*realSpaceCutoff=*/12.0, /*kMax=*/8);
    const aleator::core::ParticleData particles;
    const aleator::core::Lattice lattice;
    const EmptyNeighborList neighbors;

    REQUIRE_THROWS_AS(ewald.computeEnergy(particles, lattice, neighbors),
                       aleator::NotImplemented);
}
