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

**Status: open, unresolved, owned.** Not closed by widening any tolerance.
