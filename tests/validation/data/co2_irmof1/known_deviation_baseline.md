# Known-deviation baseline: CO2/IRMOF-1 vs. published CRAFTED reference

Per CLAUDE.md section 4. Documents where the baseline in
`tests/known_deviation/test_gcmc_co2_irmof1_known_deviation.cc` came from.

## Provenance

Recorded from a real calibration run of the exact GCMC protocol these
tests use (seed 2026, streamIndex = point index, 2,000-step equilibration
+ 40,000-step production, 8 blocks, `dev` preset, macOS arm64, 2026-08-04).

| P (Pa) | computed (mol/kg) | reference, CRAFTED (mol/kg) | combined sigma | relative deviation |
|---|---|---|---|---|
| 1e4 | 0.0753 +/- 0.0035 | 0.0829 +/- 0.0817 | 0.09 | 9.2% |
| 1e5 | 0.8767 +/- 0.0467 | 0.8304 +/- 0.2343 | 0.19 | 5.6% |
| 1e6 | 10.3293 +/- 0.5966 | 12.8323 +/- 0.9939 | **2.16** | **19.5%** |

The first two points are well within statistical noise (both under 0.2
combined sigma) and are asserted tightly in
`tests/validation/test_gcmc_co2_irmof1_isotherm.cc` (gate: 2.0 combined
sigma, >10x margin over what was actually observed). The third
(highest-pressure) point shows a real, non-trivial gap and is tracked
here instead, per section 4, rather than folded into the tight test
behind a loosened tolerance.

## What's independently verified

The electrostatic (reciprocal/self/exclusion) cache is correct at this
exact system size (848 framework atoms, 2x2x2 supercell) and occupancy
(up to ~127 molecules): every point in both the tight test and this
known-deviation test checks the incrementally-maintained Ewald state
against a full from-scratch recomputation and passes at machine precision
(relative drift 0.00e+00 to 1.79e-16 across all three calibration points).
The gap at 1e6 Pa is therefore not attributable to a bug in the
charged-GCMC energy bookkeeping (`EwaldIncrementalState`,
`MonteCarloEngine`'s charged move integration) -- see
`tests/validation/test_ewald_incremental_state.cc` and
`tests/validation/test_gcmc_charged_detailed_balance.cc` for the
dedicated tests of that mechanism.

## Leading (not independently confirmed) hypothesis

`kEquilibrationSteps = 2000` was chosen so the tight test's two points run
in a practical amount of time within this session's budget. It is
calibrated against this test's own low occupancy (<=~11 molecules at 1e4
and 1e5 Pa) and was very plausibly insufficient for the 1e6 Pa point's
much higher equilibrium occupancy (~127 molecules) to actually be reached
from an empty starting configuration within only 2000 steps -- an
under-equilibration artifact of this session's time budget, not
necessarily a force-field or sampling defect.

**Not verified.** A materially longer equilibration at ~127-molecule
occupancy was not completed in this session: this system's per-step cost
scales with occupancy (each trial move's real-space Coulomb/dispersion
scan is O(total charges present)), and a long run at this occupancy
already took on the order of tens of minutes for the production phase
alone during calibration. Confirming or ruling out this hypothesis (e.g.
by re-running with 10x-20x more equilibration steps, or by comparing
against a longer-equilibrated Widom-insertion self-consistency check
analogous to the methane/IRMOF-1 tight test) is a real, scoped, owned
follow-up -- not attempted here.

## Other hypotheses not yet checked

Unlike the methane/IRMOF-1 investigation, this session did not have
remaining time budget to check: whether CRAFTED's own RASPA run used a
materially different effective Ewald precision/kMax than this codebase's
kMax=6 (a real convergence check was done on the *static* framework
energy, see `test_gcmc_co2_irmof1_isotherm.cc`'s comments, but not
specifically on the *loading* at 1e6 Pa's high occupancy, where
convergence behavior could differ); and whether RASPA's `EwaldPrecision
1.0e-6`-derived alpha/kMax combination differs enough from this
codebase's explicit alpha=0.25/kMax=6 to matter at high occupancy
specifically. Both are real, open, scoped follow-ups.

**Status: open, unresolved, owned.** Not closed by widening any
tolerance.
