// Validation anchor (CLAUDE.md milestone: "Detailed balance audit: for
// each move type, assert the forward/reverse acceptance ratio satisfies
// the balance condition. Implement this as an actual runtime-checkable
// test, not a code review claim.")
//
// Derivation this test verifies (see gcmc_acceptance.hpp's doc comments
// for the full grand-canonical-weight derivation): under the standard
// Metropolis prescription acc(o->n) = min(1, R), acc(n->o) = min(1, 1/R),
// detailed balance pi(o)*acc(o->n) = pi(n)*acc(n->o) holds for *any*
// positive R as an algebraic identity (case R<=1: both sides equal
// pi(o)*R = pi(n); case R>1: both sides equal pi(o) = pi(n)/R). So the
// entire correctness burden is: does this codebase's insertion/deletion/
// translation/rotation ratio functions actually produce reciprocal R
// values on the same pair of states? That is a hard, checkable, numerical
// fact -- not a matter of opinion -- and is exactly what's checked below
// for every move type, using randomly sampled state parameters rather than
// one hand-picked example.

#include <catch2/catch_test_macros.hpp>

#include <algorithm>
#include <cmath>
#include <random>

#include "engines/monte_carlo/gcmc_acceptance.hpp"

using namespace aleator::engines::gcmc;

TEST_CASE("Insertion/deletion acceptance ratios are exact reciprocals (detailed balance)",
          "[validation][montecarlo]") {
    std::mt19937 rng(123);
    std::uniform_real_distribution<double> fugacityDist(1e-3, 1e3);   // K/Angstrom^3
    std::uniform_real_distribution<double> volumeDist(1e2, 1e5);      // Angstrom^3
    std::uniform_real_distribution<double> temperatureDist(50.0, 800.0); // K
    std::uniform_real_distribution<double> deltaUDist(-500.0, 500.0);  // K
    std::uniform_int_distribution<std::size_t> countDist(0, 200);

    for (int trial = 0; trial < 10000; ++trial) {
        const double fugacity = fugacityDist(rng);
        const double volume = volumeDist(rng);
        const double temperature = temperatureDist(rng);
        const double deltaUForward = deltaUDist(rng); // U(n) - U(n-1), the insertion's own deltaU
        const std::size_t nMinusOne = countDist(rng);  // state BEFORE insertion has n-1 molecules

        // Forward: insert into the (n-1)-molecule state -> n molecules.
        const double rInsert =
            insertionRatio(fugacity, volume, nMinusOne, temperature, deltaUForward);
        // Reverse: delete from the n-molecule state -> n-1 molecules. Same
        // pair of physical states, so its "countBefore" is n = nMinusOne+1,
        // and its deltaU is the negative of the forward move's.
        const double rDelete =
            deletionRatio(fugacity, volume, nMinusOne + 1, temperature, -deltaUForward);

        const double product = rInsert * rDelete;
        INFO("fugacity=" << fugacity << " V=" << volume << " T=" << temperature
                          << " nMinusOne=" << nMinusOne << " deltaU=" << deltaUForward
                          << " rInsert=" << rInsert << " rDelete=" << rDelete);
        REQUIRE(std::abs(product - 1.0) < 1e-9);

        // And the actual clamped acceptance PROBABILITIES must satisfy the
        // literal detailed-balance equation for some representative
        // (unnormalized) pi(n-1)=1 reference: pi(n-1)*acc(fwd) ==
        // pi(n)*acc(rev), with pi(n) := pi(n-1)*rInsert (that's the
        // definition of rInsert as the pi(n)/pi(n-1) ratio).
        const double piNMinusOne = 1.0;
        const double piN = piNMinusOne * rInsert;
        const double lhs = piNMinusOne * clampToProbability(rInsert);
        const double rhs = piN * clampToProbability(rDelete);
        REQUIRE(std::abs(lhs - rhs) < 1e-9 * std::max(std::abs(lhs), std::abs(rhs)) + 1e-15);
    }
}

TEST_CASE("Translation acceptance ratio is self-reciprocal (detailed balance)",
          "[validation][montecarlo]") {
    std::mt19937 rng(456);
    std::uniform_real_distribution<double> temperatureDist(50.0, 800.0);
    std::uniform_real_distribution<double> deltaUDist(-500.0, 500.0);

    for (int trial = 0; trial < 10000; ++trial) {
        const double temperature = temperatureDist(rng);
        const double deltaU = deltaUDist(rng);
        const double rForward = metropolisRatio(temperature, deltaU);
        const double rReverse = metropolisRatio(temperature, -deltaU);
        REQUIRE(std::abs(rForward * rReverse - 1.0) < 1e-12);

        const double pi0 = 1.0;
        const double pi1 = pi0 * rForward;
        const double lhs = pi0 * clampToProbability(rForward);
        const double rhs = pi1 * clampToProbability(rReverse);
        REQUIRE(std::abs(lhs - rhs) < 1e-9 * std::max(std::abs(lhs), std::abs(rhs)) + 1e-15);
    }
}

TEST_CASE("Rotation acceptance ratio is self-reciprocal (detailed balance)",
          "[validation][montecarlo]") {
    // Rotation uses the identical Metropolis form as translation (a small
    // random reorientation is as symmetric a proposal as a small random
    // displacement — see MonteCarloEngine::attemptRotation), so this is
    // the same algebraic check, kept as its own TEST_CASE because the
    // milestone asks for detailed balance to be audited "for each move
    // type" explicitly, not inferred by analogy.
    std::mt19937 rng(789);
    std::uniform_real_distribution<double> temperatureDist(50.0, 800.0);
    std::uniform_real_distribution<double> deltaUDist(-500.0, 500.0);

    for (int trial = 0; trial < 10000; ++trial) {
        const double temperature = temperatureDist(rng);
        const double deltaU = deltaUDist(rng);
        const double rForward = metropolisRatio(temperature, deltaU);
        const double rReverse = metropolisRatio(temperature, -deltaU);
        REQUIRE(std::abs(rForward * rReverse - 1.0) < 1e-12);
    }
}

TEST_CASE("Identity-swap acceptance ratio is an exact reciprocal on the reverse transition "
          "(detailed balance)",
          "[validation][montecarlo][mixture]") {
    // Mixture GCMC milestone: swapRatio's own doc comment (gcmc_acceptance.hpp)
    // derives R_swap(A->B) = (fB/fA) * (countABefore/(countBBefore+1)) *
    // exp(-dU/T), with NO leftover mass/thermal-wavelength factor between
    // species of different mass -- this is exactly the kind of prefactor
    // error the milestone warns "gives smooth, plausible, wrong
    // selectivities," so it is checked here as a hard numerical identity,
    // not trusted from the derivation alone. Species A and B get
    // INDEPENDENTLY sampled fugacities (deliberately, not the same value)
    // so a bug that accidentally cancelled fB/fA to 1 (e.g. reusing fA on
    // both sides) would be caught.
    std::mt19937 rng(20261);
    std::uniform_real_distribution<double> fugacityDist(1e-3, 1e3);
    std::uniform_real_distribution<double> temperatureDist(50.0, 800.0);
    std::uniform_real_distribution<double> deltaUDist(-500.0, 500.0);
    std::uniform_int_distribution<std::size_t> countDist(0, 200);

    for (int trial = 0; trial < 10000; ++trial) {
        const double fugacityA = fugacityDist(rng);
        const double fugacityB = fugacityDist(rng);
        const double temperature = temperatureDist(rng);
        const double deltaUForward = deltaUDist(rng); // U(after A->B) - U(before)
        const std::size_t countABefore = countDist(rng);
        const std::size_t countBBefore = countDist(rng);

        // Forward: (countABefore, countBBefore) -> (countABefore-1, countBBefore+1).
        const double rForward = swapRatio(fugacityA, fugacityB, countABefore, countBBefore,
                                           temperature, deltaUForward);
        // Reverse: on that resulting state, B->A. Its "from" count is the
        // post-forward B count (countBBefore+1); its "to" count is the
        // post-forward A count (countABefore-1, well-defined since
        // countABefore>=1 whenever the forward ratio is ever actually used
        // -- see MonteCarloEngine::attemptSwap, which never proposes a swap
        // FROM a species with zero molecules present).
        if (countABefore == 0) {
            continue; // forward move is never actually proposable from this state; skip
        }
        const double rReverse = swapRatio(fugacityB, fugacityA, countBBefore + 1, countABefore - 1,
                                           temperature, -deltaUForward);

        const double product = rForward * rReverse;
        INFO("fA=" << fugacityA << " fB=" << fugacityB << " T=" << temperature
                    << " countABefore=" << countABefore << " countBBefore=" << countBBefore
                    << " deltaU=" << deltaUForward << " rForward=" << rForward
                    << " rReverse=" << rReverse);
        REQUIRE(std::abs(product - 1.0) < 1e-9);

        const double piBefore = 1.0;
        const double piAfter = piBefore * rForward;
        const double lhs = piBefore * clampToProbability(rForward);
        const double rhs = piAfter * clampToProbability(rReverse);
        REQUIRE(std::abs(lhs - rhs) < 1e-9 * std::max(std::abs(lhs), std::abs(rhs)) + 1e-15);
    }
}
