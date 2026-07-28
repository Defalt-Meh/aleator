# Source

## Crystal structure: `IRMOF-1.cif`

Fetched 2026-07-28 from the RASPA2 reference implementation's bundled MOF
structures:

- `https://raw.githubusercontent.com/numat/RASPA2/master/structures/mofs/cif/IRMOF-1.cif`
- Authored by David Dubbeldam (RASPA's original author).
- Cites the original synthesis/structure paper: Eddaoudi, M.; Kim, J.;
  Rosi, N.; Vodak, D.; Wachter, J.; O'Keeffe, M.; Yaghi, O. M. "Systematic
  design of pore size and functionality in isoreticular MOFs and their
  application in methane storage." *Science* **2002**, *295* (5554),
  469-472.
- Cubic, space group Fm-3m (#225), a = 25.832 Angstrom, 7 asymmetric-unit
  atoms (Zn, 2xO, 3xC, H) expanded via this codebase's own CIF symmetry
  machinery (io/, validated in an earlier milestone against real IZA
  zeolite structures) — not pre-expanded or hand-edited.
- All `_atom_site_charge` values in the file are 0 — RASPA's own
  "GenericMOFs" force field (see below) treats IRMOF-1 as a non-polar,
  LJ-only framework, so this codebase's GCMC validation does not need
  electrostatics (Ewald) wired in for this system.

## Force field parameters (framework + methane)

Fetched 2026-07-28 from RASPA2's bundled "GenericMOFs" (UFF-derived) and
species definitions:

- `https://raw.githubusercontent.com/numat/RASPA2/master/forcefield/GenericMOFs/pseudo_atoms.def`
- `https://raw.githubusercontent.com/numat/RASPA2/master/forcefield/GenericMOFs/force_field_mixing_rules.def`

Transcribed values used directly (epsilon in K, sigma in Angstrom, already
in this codebase's internal unit convention — no conversion needed):

| Species | epsilon/k_B (K) | sigma (Angstrom) | mass (amu) |
|---|---|---|---|
| Zn (framework) | 62.3992 | 2.46155 | 65.38  |
| O  (framework) | 48.1581 | 3.03315 | 16.00  |
| C  (framework) | 47.8562 | 3.47299 | 12.01  |
| H  (framework) |  7.64893 | 2.84642 | 1.008 |
| CH4_sp3 (adsorbate, united-atom) | 158.5 | 3.72 | 16.04246 |

Mixing rule: Lorentz-Berthelot (arithmetic mean sigma, geometric mean
epsilon) — RASPA's file states this explicitly ("general mixing rule for
Lennard-Jones: Lorentz-Berthelot"), matching this codebase's existing
`MixingRule::LorentzBerthelot` (validated in the forcefield/pairwise
milestone).

Truncation: RASPA's file states "shifted" with "no tail corrections" for
this force field — matches `LennardJonesTruncation::Shifted` (no
`tailEnergyCorrection` needed).

## Reference isotherm: `methane_irmof1_isotherm_298K.csv`

Fetched 2026-07-28 from the pyIAST package's bundled test data:

- `https://raw.githubusercontent.com/CorySimon/pyIAST/master/test/IRMOF-1_methane_isotherm_298K.csv`
- pyIAST: Simon, C. M.; Smit, B.; Haranczyk, M. "pyIAST: Ideal adsorbed
  solution theory (IAST) Python package." *Computer Physics
  Communications* **2016**, *200*, 364-380.
- The package's own test notebook (`test/Methane and ethane
  test.ipynb`) describes this file as "Simulated pure-component
  adsorption isotherms at 298 K" for methane in IRMOF-1 — i.e. GCMC
  simulation output (not experimental data), which is the correct kind of
  reference for validating this codebase's own GCMC implementation
  against (an experimental isotherm would also carry force-field
  transferability error on top of any simulation error, muddying which
  discrepancies are attributable to what).
- Columns: `Pressure(bar)`, `Loading(mmol/g)`. Spans 1e-4 to 150 bar.

Units conversion used when comparing against this codebase's output:
loading in molecules/Angstrom^3 -> mmol/g via the IRMOF-1 unit cell's own
mass and volume (computed from the structure this codebase's CIF reader
itself parses, not hand-entered).
