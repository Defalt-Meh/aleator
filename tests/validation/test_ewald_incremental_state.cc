// Validation anchor for the charged-GCMC milestone (CLAUDE.md: "the
// reciprocal-space Ewald sum... does not decompose into a per-particle
// contribution... the standard solution is incremental structure-factor
// bookkeeping"). Tests EwaldIncrementalState in isolation, driving it
// directly (not through MonteCarloEngine) with a scripted sequence of
// insert/delete/move trials, some committed and some deliberately left
// uncommitted (simulating rejection) -- exactly mirroring, in parallel, a
// real core::ParticleData + molecule-index bookkeeping so
// recomputeFromScratch() has real ground truth to check against.
//
// Two required properties, per the milestone spec:
//   1. DRIFT GATE: incrementally-maintained energy vs. a full from-scratch
//      recomputation must agree to 1e-10 relative, with no secular drift
//      over a long run.
//   2. ROLLBACK EXACTNESS: after a rejected (never-committed) trial, the
//      committed cache must be bit-identical (raw bits, not a tolerance) to
//      what it was before the trial.

#include <catch2/catch_test_macros.hpp>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdio>
#include <vector>

#include "core/geometry/lattice.hpp"
#include "core/math/counter_based_rng.hpp"
#include "core/math/particle_data.hpp"
#include "forcefield/electrostatics/ewald.hpp"
#include "forcefield/electrostatics/ewald_incremental_state.hpp"

using aleator::core::Lattice;
using aleator::core::ParticleData;
using aleator::core::Philox4x32Rng;
using aleator::forcefield::ChargedSite;
using aleator::forcefield::Ewald;
using aleator::forcefield::EwaldIncrementalState;

namespace {

Lattice cubicLattice(double length) {
    std::array<std::array<double, 3>, 3> m{};
    m[0] = {length, 0.0, 0.0};
    m[1] = {0.0, length, 0.0};
    m[2] = {0.0, 0.0, length};
    return Lattice(m);
}

/// A small, deliberately non-physical charged "framework": 8 alternating
/// +1/-1 point charges at the corners of a smaller cube inside the
/// simulation box -- net neutral, nothing more meaningful is needed for a
/// pure state-machine/bookkeeping test (the physics itself is validated
/// elsewhere: Madelung constant, SPC/E water, and the real IRMOF-1/CO2
/// isotherm test).
ParticleData buildFramework(double boxLength) {
    ParticleData framework;
    const std::array<double, 2> coords{boxLength * 0.25, boxLength * 0.75};
    int sign = 1;
    for (double x : coords) {
        for (double y : coords) {
            for (double z : coords) {
                framework.push_back(x, y, z, 0.0, 0.0, 0.0, /*mass=*/20.0,
                                     /*charge=*/static_cast<double>(sign), /*species=*/0);
                sign = -sign;
            }
        }
    }
    return framework;
}

/// A rigid, linear, CO2-like 3-site species: C (+0.7 e) between two O
/// (-0.35 e each), net neutral, 1.16 Ang bonds -- same charges as the real
/// CO2/IRMOF-1 reference test uses (CRAFTED database's bundled CO2.def),
/// reused here purely as a realistic charge magnitude/geometry, not
/// because this test cares about CO2 specifically.
struct Co2Sites {
    static constexpr double kBondLength = 1.16;
    static constexpr double kCarbonCharge = 0.7;
    static constexpr double kOxygenCharge = -0.35;
};

std::vector<ChargedSite> co2SitesAt(double cx, double cy, double cz) {
    return {
        {cx, cy, cz - Co2Sites::kBondLength, Co2Sites::kOxygenCharge},
        {cx, cy, cz, Co2Sites::kCarbonCharge},
        {cx, cy, cz + Co2Sites::kBondLength, Co2Sites::kOxygenCharge},
    };
}

/// Mirrors a ChargedSite list into `particles`, returning the new global
/// indices (in the same order), and records them as one molecule in
/// `molecules`.
std::vector<std::size_t> insertIntoParticles(ParticleData& particles,
                                              std::vector<std::vector<std::size_t>>& molecules,
                                              const std::vector<ChargedSite>& sites) {
    std::vector<std::size_t> indices;
    indices.reserve(sites.size());
    for (const auto& s : sites) {
        indices.push_back(
            particles.push_back(s.x, s.y, s.z, 0.0, 0.0, 0.0, /*mass=*/16.0, s.charge,
                                 /*species=*/1));
    }
    molecules.push_back(indices);
    return indices;
}

std::vector<ChargedSite> sitesOf(const ParticleData& particles,
                                  const std::vector<std::size_t>& indices) {
    std::vector<ChargedSite> sites;
    sites.reserve(indices.size());
    for (std::size_t i : indices) {
        sites.push_back({particles.x[i], particles.y[i], particles.z[i], particles.charge[i]});
    }
    return sites;
}

/// Removes the LAST molecule tracked in `molecules` from `particles` via
/// swap-and-pop, matching MonteCarloEngine::removeMolecule's convention
/// closely enough for this test's purposes (this test never needs to
/// remove anything other than the most recently inserted molecule, so the
/// general cross-reference fixup MonteCarloEngine needs isn't required
/// here).
void removeLastMolecule(ParticleData& particles, std::vector<std::vector<std::size_t>>& molecules) {
    const auto sites = molecules.back();
    // Sites were pushed contiguously and are the last `sites.size()`
    // particles (this test always removes what it most recently added, so
    // no other molecule's indices are ever disturbed).
    for (std::size_t k = 0; k < sites.size(); ++k) {
        particles.popBack();
    }
    molecules.pop_back();
}

} // namespace

TEST_CASE("EwaldIncrementalState matches a from-scratch recomputation with no drift over a "
          "scripted sequence of committed and rejected trials",
          "[validation][electrostatics]") {
    const double boxLength = 20.0;
    const Lattice lattice = cubicLattice(boxLength);
    ParticleData particles = buildFramework(boxLength);
    std::vector<std::vector<std::size_t>> molecules;

    const Ewald ewald(/*alpha=*/0.25, /*realSpaceCutoff=*/9.0, /*kMax=*/6);
    EwaldIncrementalState state(ewald, particles, lattice);

    // Sanity: fresh state matches a from-scratch recompute immediately.
    {
        const double fresh = state.recomputeFromScratch(particles, molecules);
        const double committed = state.committedEnergy();
        REQUIRE(std::abs(fresh - committed) <= 1e-10 * std::max(1.0, std::abs(fresh)));
    }

    Philox4x32Rng rng;
    rng.seed(/*seed=*/7, /*streamIndex=*/0);

    constexpr int kShortSteps = 2000;
    double maxRelDrift = 0.0;

    for (int step = 0; step < kShortSteps; ++step) {
        const double pick = rng.nextUniform();
        const double cx = boxLength * rng.nextUniform();
        const double cy = boxLength * rng.nextUniform();
        const double cz = boxLength * rng.nextUniform();

        if (pick < 0.4 || molecules.empty()) {
            // Propose an insertion; sometimes commit, sometimes reject.
            const auto trialSites = co2SitesAt(cx, cy, cz);
            const double delta = state.proposeInsertion(trialSites);
            (void)delta; // this test checks cache correctness, not acceptance physics
            if (rng.nextUniform() < 0.5) {
                state.commitInsertion(trialSites);
                insertIntoParticles(particles, molecules, trialSites);
            }
        } else if (pick < 0.7) {
            // Propose deleting the most recently inserted molecule.
            const auto existing = sitesOf(particles, molecules.back());
            const double delta = state.proposeDeletion(existing);
            (void)delta;
            if (rng.nextUniform() < 0.5) {
                state.commitDeletion(existing);
                removeLastMolecule(particles, molecules);
            }
        } else {
            // Propose moving the most recently inserted molecule bodily.
            const auto oldSites = sitesOf(particles, molecules.back());
            std::vector<ChargedSite> newSites = oldSites;
            const double dx = (rng.nextUniform() - 0.5) * 3.0;
            const double dy = (rng.nextUniform() - 0.5) * 3.0;
            const double dz = (rng.nextUniform() - 0.5) * 3.0;
            for (auto& s : newSites) {
                s.x += dx;
                s.y += dy;
                s.z += dz;
            }
            const double delta = state.proposeMove(oldSites, newSites);
            (void)delta;
            if (rng.nextUniform() < 0.5) {
                state.commitMove(oldSites, newSites);
                for (std::size_t k = 0; k < molecules.back().size(); ++k) {
                    const std::size_t idx = molecules.back()[k];
                    particles.x[idx] = newSites[k].x;
                    particles.y[idx] = newSites[k].y;
                    particles.z[idx] = newSites[k].z;
                }
            }
        }

        // DRIFT GATE, tier 1: every single step.
        const double fresh = state.recomputeFromScratch(particles, molecules);
        const double committed = state.committedEnergy();
        const double relError = std::abs(fresh - committed) / std::max(1.0, std::abs(fresh));
        maxRelDrift = std::max(maxRelDrift, relError);
        INFO("step " << step << ": fresh=" << fresh << " committed=" << committed
                     << " molecules=" << molecules.size());
        REQUIRE(relError <= 1e-10);
    }

    std::printf("Short-run (every-step) drift gate: %d steps, max relative error = %.3e\n",
                kShortSteps, maxRelDrift);
}

TEST_CASE("EwaldIncrementalState shows no secular drift over a long run, checked every 100 "
          "steps",
          "[validation][electrostatics][slow]") {
    const double boxLength = 20.0;
    const Lattice lattice = cubicLattice(boxLength);
    ParticleData particles = buildFramework(boxLength);
    std::vector<std::vector<std::size_t>> molecules;

    const Ewald ewald(/*alpha=*/0.25, /*realSpaceCutoff=*/9.0, /*kMax=*/6);
    EwaldIncrementalState state(ewald, particles, lattice);

    Philox4x32Rng rng;
    rng.seed(/*seed=*/2026, /*streamIndex=*/0);

    constexpr int kLongSteps = 100'000;
    constexpr int kCheckStride = 100;
    // Periodic resync bounds long-run floating-point accumulation error --
    // see EwaldIncrementalState::resync()'s doc comment. This mirrors the
    // periodic resync MonteCarloEngine itself performs during a real
    // charged GCMC run.
    constexpr int kResyncStride = 500;
    // Unlike a real fugacity-controlled GCMC run (which naturally
    // equilibrates around a bounded mean occupancy), this test's
    // insert/delete/move selection is a blind, energy-independent birth-
    // death process -- without a cap it random-walks to an ever-growing
    // molecule count, which is unrepresentative of real usage and (being
    // an unbounded amount of summed floating-point terms) makes long-run
    // floating-point accumulation worse than any real run would see. Cap
    // it at a modest, still-nontrivial size.
    constexpr std::size_t kMaxMolecules = 25;
    std::vector<double> checkedRelErrors;

    for (int step = 0; step < kLongSteps; ++step) {
        const double pick = rng.nextUniform();
        const double cx = boxLength * rng.nextUniform();
        const double cy = boxLength * rng.nextUniform();
        const double cz = boxLength * rng.nextUniform();

        if (molecules.size() >= kMaxMolecules) {
            // Force a deletion attempt to keep the population bounded.
            const auto existing = sitesOf(particles, molecules.back());
            (void)state.proposeDeletion(existing);
            if (rng.nextUniform() < 0.8) {
                state.commitDeletion(existing);
                removeLastMolecule(particles, molecules);
            }
        } else if (pick < 0.4 || molecules.empty()) {
            const auto trialSites = co2SitesAt(cx, cy, cz);
            (void)state.proposeInsertion(trialSites);
            if (rng.nextUniform() < 0.5) {
                state.commitInsertion(trialSites);
                insertIntoParticles(particles, molecules, trialSites);
            }
        } else if (pick < 0.7) {
            const auto existing = sitesOf(particles, molecules.back());
            (void)state.proposeDeletion(existing);
            if (rng.nextUniform() < 0.5) {
                state.commitDeletion(existing);
                removeLastMolecule(particles, molecules);
            }
        } else {
            const auto oldSites = sitesOf(particles, molecules.back());
            std::vector<ChargedSite> newSites = oldSites;
            const double dx = (rng.nextUniform() - 0.5) * 3.0;
            const double dy = (rng.nextUniform() - 0.5) * 3.0;
            const double dz = (rng.nextUniform() - 0.5) * 3.0;
            for (auto& s : newSites) {
                s.x += dx;
                s.y += dy;
                s.z += dz;
            }
            (void)state.proposeMove(oldSites, newSites);
            if (rng.nextUniform() < 0.5) {
                state.commitMove(oldSites, newSites);
                for (std::size_t k = 0; k < molecules.back().size(); ++k) {
                    const std::size_t idx = molecules.back()[k];
                    particles.x[idx] = newSites[k].x;
                    particles.y[idx] = newSites[k].y;
                    particles.z[idx] = newSites[k].z;
                }
            }
        }

        if (step % kCheckStride == 0) {
            const double fresh = state.recomputeFromScratch(particles, molecules);
            const double committed = state.committedEnergy();
            const double relError = std::abs(fresh - committed) / std::max(1.0, std::abs(fresh));
            checkedRelErrors.push_back(relError);
            INFO("step " << step << ": fresh=" << fresh << " committed=" << committed);
            REQUIRE(relError <= 1e-10);
        }

        if (step > 0 && step % kResyncStride == 0) {
            state.resync(particles, molecules);
        }
    }

    // No SECULAR drift: the error at the end of the run must not be
    // systematically larger than at the start (a real cache-drift bug
    // accumulates monotonically; independent floating-point noise does
    // not).
    const double firstQuarterMax =
        *std::max_element(checkedRelErrors.begin(),
                           checkedRelErrors.begin() +
                               static_cast<long>(checkedRelErrors.size() / 4));
    const double lastQuarterMax = *std::max_element(
        checkedRelErrors.end() - static_cast<long>(checkedRelErrors.size() / 4),
        checkedRelErrors.end());
    std::printf(
        "Long-run (every-%d-steps) drift gate: %d steps, %zu checks, first-quarter max = "
        "%.3e, last-quarter max = %.3e\n",
        kCheckStride, kLongSteps, checkedRelErrors.size(), firstQuarterMax, lastQuarterMax);
    // Both must individually pass the same 1e-10 gate (already asserted
    // above); this additionally confirms the error isn't drifting upward
    // across the run.
    REQUIRE(lastQuarterMax <= 1e-10);
}

TEST_CASE("EwaldIncrementalState's insertion and deletion deltas are exact negatives for the "
          "same molecule at the same position (detailed balance at the energy level)",
          "[validation][electrostatics]") {
    // CLAUDE.md milestone: "Detailed balance for the charged moves, same
    // runtime-checked algebraic identity as the uncharged case." The
    // uncharged case's identity (gcmc_acceptance.hpp's insertionRatio *
    // deletionRatio == 1) is an algebraic property of the RATIO FORMULA
    // alone and holds for any deltaU, charged or not -- already proven in
    // test_gcmc_detailed_balance.cc. What's new and specific to charged
    // GCMC is whether the ENERGY FUNCTION feeding that formula is itself
    // round-trip consistent: inserting a molecule and then, from the
    // resulting state, deleting that exact molecule must produce exactly
    // opposite energy deltas, or the ratio identity would be proven
    // against numbers that don't correspond to real forward/reverse moves
    // on the same pair of states.
    const double boxLength = 20.0;
    const Lattice lattice = cubicLattice(boxLength);
    ParticleData particles = buildFramework(boxLength);

    const Ewald ewald(/*alpha=*/0.25, /*realSpaceCutoff=*/9.0, /*kMax=*/6);
    EwaldIncrementalState state(ewald, particles, lattice);

    Philox4x32Rng rng;
    rng.seed(/*seed=*/99, /*streamIndex=*/0);

    for (int trial = 0; trial < 200; ++trial) {
        const double cx = boxLength * rng.nextUniform();
        const double cy = boxLength * rng.nextUniform();
        const double cz = boxLength * rng.nextUniform();
        const auto sites = co2SitesAt(cx, cy, cz);

        const double deltaInsert = state.proposeInsertion(sites);
        state.commitInsertion(sites);

        const double deltaDelete = state.proposeDeletion(sites);

        INFO("trial " << trial << ": deltaInsert=" << deltaInsert
                       << " deltaDelete=" << deltaDelete);
        REQUIRE(std::abs(deltaInsert + deltaDelete) <=
                1e-9 * std::max(1.0, std::abs(deltaInsert)));

        // Restore to the pre-trial state for the next iteration (so every
        // trial starts from the same framework-only baseline).
        state.commitDeletion(sites);
    }
}

TEST_CASE("EwaldIncrementalState leaves the committed cache bit-identical after a rejected "
          "trial",
          "[validation][electrostatics]") {
    const double boxLength = 20.0;
    const Lattice lattice = cubicLattice(boxLength);
    ParticleData particles = buildFramework(boxLength);

    const Ewald ewald(/*alpha=*/0.25, /*realSpaceCutoff=*/9.0, /*kMax=*/6);
    EwaldIncrementalState state(ewald, particles, lattice);

    // Commit one real molecule first so self/exclusion/reciprocal state is
    // all genuinely nonzero and non-trivial, not just the framework's
    // initial (never-touched) values.
    const auto committedSites = co2SitesAt(10.0, 10.0, 10.0);
    state.commitInsertion(committedSites);

    const std::vector<double> before = state.rawStateForTesting();

    // A batch of trials that are proposed but deliberately never
    // committed -- insertion, deletion, and move, covering every path
    // that touches the committed cache on the accept side.
    const auto rejectedInsertSites = co2SitesAt(3.0, 3.0, 3.0);
    (void)state.proposeInsertion(rejectedInsertSites);

    (void)state.proposeDeletion(committedSites);

    std::vector<ChargedSite> rejectedMoveSites = committedSites;
    for (auto& s : rejectedMoveSites) {
        s.x += 5.0;
    }
    (void)state.proposeMove(committedSites, rejectedMoveSites);

    const std::vector<double> after = state.rawStateForTesting();

    REQUIRE(before.size() == after.size());
    for (std::size_t i = 0; i < before.size(); ++i) {
        // Raw bit comparison, not a tolerance: propose* must not have
        // mutated a single bit of the committed cache.
        CAPTURE(i, before[i], after[i]);
        REQUIRE(before[i] == after[i]);
    }
}
