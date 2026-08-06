# Known-deviation baseline: CH4/IRMOF-1 vs. published pyIAST isotherm

Per CLAUDE.md section 4 ("Handling a validated component that disagrees
with its reference"). This file documents where the baseline numbers in
`test_gcmc_ch4_irmof1_isotherm.cc`'s known-deviation test case came from —
it is not parsed at runtime (the numbers are transcribed into the test,
same convention as `methane_irmof1_isotherm_298K.csv` and
`referenceIsotherm()`).

## Provenance

Recorded from a real run of this exact test (seed 2026, streamIndex =
point index, 30,000-step equilibration + 90,000-step production, 10
blocks) under the `asan` preset (UBSan; ASan disabled) on 2026-08-02. Per
CLAUDE.md invariant #5 (deterministic, reproducible RNG), this run's
numbers were independently confirmed bit-identical to an earlier `dev`
preset run of the same test on the same day.

| P (bar) | computed (mmol/g) | reference, pyIAST (mmol/g) | relative deviation |
|---|---|---|---|
| 0.1  | 0.0376 | 0.04301247025401385 | 12.58% |
| 1.0  | 0.3860 | 0.4384412654666955  | 11.96% |
| 5.0  | 1.9061 | 2.2082018401014736  | 13.68% |
| 10.0 | 3.7551 | 4.389245423072791   | 14.45% |

## Sampling machinery is not the explanation

`tests/validation/test_gcmc_ch4_irmof1_isotherm.cc` (the tight, split-off
self-consistency test) measured this engine's own low-pressure loading on
the real IRMOF-1 structure against its own Widom-insertion Henry
coefficient on the same structure -- no external reference involved at
all. Result: relative difference 3.42%, block-averaged relative standard
error 1.49% (~2.3 standard errors, consistent with ordinary Monte Carlo
noise). That is far tighter than the 12-15% gap against the external
pyIAST reference, which is itself evidence that the gap is a genuine
external/force-field/reference discrepancy rather than a bug in the
acceptance criteria, the force-field energy call, or the fugacity
conversion.

## Investigation performed this session (CLAUDE.md hygiene milestone, section 0 defect 3)

Hypotheses checked (see CLAUDE.md section 0 and this test file's own
comments for the full writeup):

- **Tail-correction convention**: RULED OUT as the explanation. RASPA2's
  GenericMOFs force-field file (`force_field_mixing_rules.def`, fetched
  live and re-verified during this session) states "shifted" truncation
  and "no tailcorrections" explicitly — exactly what this codebase already
  does (`LennardJonesTruncation::Shifted`, no `tailEnergyCorrection` term
  added). The standing hypothesis in earlier session notes was that a
  *missing* tail correction explained the shortfall; the force field's own
  stated convention says there should be no tail correction to add in the
  first place, so this is not an available explanation on its own.
- **Mixing rule / LJ parameter set**: CONFIRMED correct. Live-refetched
  `pseudo_atoms.def` and `force_field_mixing_rules.def` from RASPA2 (not
  re-trusted from an earlier session's cached transcription) match this
  test's `uffParameters()` table exactly (Zn, O, C, H epsilon/sigma) and
  confirm Lorentz-Berthelot mixing, matching `MixingRule::LorentzBerthelot`
  (this codebase's default, used here).
- **Cutoff / minimum-image margin**: RULED OUT, for the within-cell portion
  of this hypothesis. The 12.0 Angstrom cutoff leaves only ~0.9 Angstrom
  of margin under minimum image in this 25.832 Angstrom cubic cell (L/2 =
  12.916 Ang) — safe (no double-counting), but tight. A real single-point
  sensitivity check was run (P=1.0 bar, cutoff 12.0 vs. 12.8 Ang, same
  unit cell, same seed/statistics protocol as the main test): cutoff=12.0
  gave `<N> = 2.4305 +/- 0.0649`; cutoff=12.8 gave `<N> = 2.4470 +/-
  0.0645` — a difference of 0.25 standard errors, i.e. statistically
  indistinguishable. The dispersion energy in the last ~0.8 Angstrom
  shell before the minimum-image limit does not contribute meaningfully,
  consistent with 1/r^6 falloff. This does NOT rule out a genuinely larger
  cutoff enabled by a bigger simulation cell (a single 25.832 Angstrom
  cubic cell has no additional periodic images to capture beyond L/2
  regardless of the nominal cutoff value) — a full 2x2x2 supercell with a
  substantially larger cutoff (which would test whether pyIAST's
  underlying RASPA run used a bigger simulation cell) was NOT attempted:
  this codebase's GCMC engine bypasses `core/neighbor` in favor of an O(N)
  scan (CLAUDE.md section 0, defect 4), and an 8x larger framework (~3,392
  atoms) would make that prohibitively slow for this session's time
  budget. Left as a genuinely open, scoped follow-up, not resolved.
- **pyIAST's own simulation parameters**: pyIAST's public repository
  (`CorySimon/pyIAST`, `test/` directory and its two notebooks) documents
  the reference CSV as "simulated pure-component adsorption isotherms" but
  does NOT state the cutoff, supercell size, or exact force field used to
  generate them — there is no further public documentation to confirm or
  rule out a simulation-setup mismatch against. This is itself a finding:
  part of the gap is unverifiable from available sources, not just
  unresolved.
- **Inaccessible pore blocking**: RULED OUT by reading this codebase's own
  insertion code (`MonteCarloEngine::attemptInsertion`,
  `src/engines/monte_carlo/monte_carlo_engine.cc`): trial insertions are
  sampled uniformly over the full unit cell in fractional coordinates,
  with no pore-accessibility filtering at all. An insertion attempt inside
  a genuinely inaccessible pocket is simply rejected by the ordinary
  Metropolis criterion (near-zero acceptance probability from severe
  overlap energy) — the same outcome a reference simulation gets by
  excluding that volume from attempts in the first place. Blocking or not
  blocking changes sampling efficiency, not the equilibrium loading, so it
  cannot explain a systematic several-percent shortfall.

## Investigation performed this session (supercell replication milestone)

The prior session left a real cutoff/supercell test genuinely un-attempted
("prohibitively slow" for the O(N) guest-host scan at the time). This
session's own milestone (automatic supercell replication,
`core::minimumSupercellReplication`/`core::replicateSupercell`, see CLAUDE.md
section 0) made it possible for the first time and it was actually run, on
a Release build against real IRMOF-1 and this exact test's own 0.1 bar
point (same seed=2026, same 30,000+90,000 step protocol).

**Accumulated shifted-potential offset, computed directly.** The standing
hypothesis was: a 13% near-constant ratio across two decades of pressure
looks like an energy-scale offset of ~36 K per molecule
(`ln(1.13)*298 K`), the size of an accumulated `-V_LJ(r_c)` (shifted-LJ
truncation constant) term. Computed directly rather than estimated: the
real CH4-framework LJ energy was minimized over the real IRMOF-1 unit cell
(30^3 coarse grid + local refinement) to find the actual lowest-energy
adsorption site (E = -1696.3 K at fractional (0.150, 0.153, 0.850)), then
`sum(-V_LJ(cutoff=12.0))` was computed over every one of the 176 framework
atoms within that cutoff of the site. Result: **31.33 K — 86% of the
predicted 36.42 K.** Checked for robustness, not taken from one lucky
point: the same site (by IRMOF-1's own `Fm-3m` symmetry, 192 operations)
recurs at 8 independent locations found via the same search, all giving
the identical 31.33 K. This is real, substantial, quantitative support for
the hypothesis — not proof it's the *entire* gap, but strong evidence it's
a major, real contributor.

**Direct test: does a larger cutoff (now expressible via a 2x2x2
supercell) move the computed loading?** Yes, substantially, but *not* by
simple monotonic convergence onto the pyIAST value:

| cutoff (Ang) | cell | `<N>` | loading (mmol/g) | gap to pyIAST |
|---|---|---|---|---|
| 12.0 | 1x1x1 (424 atoms) | 0.2315 +/- 0.0036 | 0.03758 +/- 0.00058 | 12.63% |
| 16.0 | 2x2x2 (3392 atoms) | 2.1448 +/- 0.0437 | 0.04353 +/- 0.00089 | **1.21%** |
| 20.0 | 2x2x2 (3392 atoms) | 2.2508 +/- 0.0384 | 0.04568 +/- 0.00078 | 6.21% |

Loading increases monotonically with cutoff (0.0376 -> 0.0435 -> 0.0457) —
a real, substantial, cutoff-driven effect, exactly the direction and
rough scale the accumulated-shift hypothesis predicts (removing more of
the artificial positive shift as the cutoff grows makes binding more
favorable, raising loading). But it does not plateau at the pyIAST value:
it passes *through* it near cutoff=16 (gap statistically close to zero:
1.21% is about 0.7 standard errors) and continues rising past it at
cutoff=20 (gap grows again, ~3.6 standard errors — a real, not
noise-level, move in the "wrong" direction relative to pyIAST). The
cutoff=16-vs-20 difference (0.00215 mmol/g) is ~2.4x the larger of the two
points' standard errors, so this is a real trend, not two noisy draws.

**Honest verdict (CLAUDE.md section 4: "if it survives, it stays owned"):**
cutoff truncation is now demonstrated, not just hypothesized, to be a real
and substantial contributor to this system's computed loading (a ~20%
swing in loading between cutoff=12 and cutoff=20) — this genuinely
explains a large fraction of why the single-cell, cutoff=12 result sits
low. But "increase the cutoff" does not cleanly converge onto pyIAST's
number; the computed curve passes near it and continues moving away,
meaning either (a) this codebase's true infinite-cutoff-limit loading is
close to but not identical to pyIAST's value at this alpha/force-field
combination, or (b) pyIAST's own (undocumented, per above) simulation
setup used a cutoff/supercell in the vicinity of what produced the
cutoff=16 near-match, which would be coincidence-shaped rather than
demonstrated. Only one pressure point (0.1 bar) and three cutoffs were
tested this session (time-scoped, not exhaustive) — the full four-point
curve at a larger cutoff, and denser cutoff sampling to characterize the
true asymptote, remain open follow-up work. **This defect is NOT closed.**
It is better characterized, with two new, real, quantitative findings
(both supporting a major real contribution from cutoff truncation) added
to the record, not swept away.

**Status: open, unresolved, owned.** Not closed by widening any tolerance.
