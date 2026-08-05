# Source

## Structure + reference isotherm: CRAFTED database

Fetched 2026-08-03 from Zenodo:

- Record: <https://zenodo.org/records/7689919>
- Archive: `CRAFTED-1.1.1.tar.xz`, DOI [10.5281/zenodo.7689919](https://doi.org/10.5281/zenodo.7689919)
- License: CC-BY-4.0
- Citation: the CRAFTED database, *"CRAFTED: An exploratory database of
  simulated adsorption isotherms of metal-organic frameworks,"*
  *Scientific Data* (2023) — real RASPA GCMC simulations of CO2/N2
  adsorption on 690 CoRE MOF 2014 structures, crossed with {UFF, DREIDING}
  force fields and six partial-charge schemes, at 273/298/323 K.

IRMOF-1 is one of the 690 structures (kept under its actual name in the
archive, `IRMOF-1.cif`, not an opaque CSD refcode). Cu-BTC/HKUST-1 (this
milestone's other suggested reference system) was checked first and
confirmed absent from this 690-structure subset (no CIF file in any of
the six charge-scheme folders contains a "Cu" atom site).

### `IRMOF-1_primitive_DDEC.cif`

Copied verbatim from `CIF_FILES/DDEC/IRMOF-1.cif`. DDEC (density-derived
electrostatic and chemical) partial charges — the most physically-grounded
of CRAFTED's six charge schemes (DFT-electron-density partitioning of the
real crystal, not a predicted/heuristic scheme). This is the **primitive**
rhombohedral cell of IRMOF-1's Fm-3m lattice (a=b=c=18.266 Angstrom,
alpha=beta=gamma=60 degrees, 106 atoms — 1/4 the volume and atom count of
the conventional cubic cell, a=25.832 Angstrom, 424 atoms, used by the
uncharged methane/IRMOF-1 validation test), computed via this codebase's
own `L_perp/2` formula to be 7.46 Angstrom — too small to support the
reference simulation's 12.8 Angstrom cutoff on its own, which is exactly
why the reference input file specifies a 2x2x2 supercell (see below).

As tabulated, the 106 atoms' charges sum to 8x10^-6 e, not exactly zero —
an ordinary DDEC-partitioning rounding residual, not a defect in this
data. `test_gcmc_co2_irmof1_isotherm.cc` redistributes this residual
uniformly across all atoms before use (standard practice for
externally-sourced partial charges; the correction is ~7.5x10^-8 e per
atom, utterly negligible next to any individual atom's charge, which is
>=1.16x10^-1 e for every species present).

### `co2_irmof1_isotherm_298K.csv`

Copied verbatim from `ISOTHERM_FILES/DDEC_IRMOF-1_UFF_CO2_298.csv`.
Columns: pressure (Pa), mean loading (mol/kg), mean error (mol/kg). 11
points, 100 Pa to 1 MPa.

### Simulation parameters (from `INPUT_FILES/DDEC_IRMOF-1_UFF_CO2_298.input`)

- `CutOffVDW` = `CutOffChargeCharge` = 12.8 Angstrom
- `ChargeMethod Ewald`, `EwaldPrecision 1.0e-6` (RASPA derives alpha/kmax
  from a precision target rather than taking them directly; this
  codebase's `Ewald` takes explicit alpha/kmax instead — see the test file
  for the values used and why they're expected to be equivalently
  converged, not identical numbers)
- `UnitCells 2 2 2` — a real 2x2x2 supercell of the primitive cell (848
  atoms, `L_perp/2` = 14.91 Angstrom), reproduced by this codebase's own
  supercell-replication helper in the test file, not approximated by
  reusing the smaller conventional cell from the methane milestone.
- Force field: `FORCEFIELDS/UFF/` bundled in the same archive —
  **bare/generic UFF** Lennard-Jones parameters (Zn: 62.38 K / 2.462
  Angstrom, C: 52.8 / 3.431, O: 30.2 / 3.118, H: 22.14 / 2.571),
  Lorentz-Berthelot mixing, "shifted" truncation with "no tail
  corrections" stated explicitly in `force_field_mixing_rules.def` — all
  **different from** RASPA2's "GenericMOFs" set used by the uncharged
  methane/IRMOF-1 test (most notably H: 22.14 K here vs. 7.65 K there;
  GenericMOFs is a tuned/modified UFF variant, this is closer to the
  literal Rappe et al. 1992 UFF values). Not interchangeable with the
  methane test's parameters — this test uses CRAFTED's own bundled set
  throughout, not a mix of the two.
- CO2: `FORCEFIELDS/UFF/CO2.def`, rigid, linear O-C-O, bond length 1.16
  Angstrom, charges C: +0.700 e, O: -0.350 e each (net neutral) — **not**
  the EPM2-style +0.6512/-0.3256 charges RASPA2's GenericMOFs bundles for
  CH4's methane comparison; a different, simpler round-number charge
  model specific to this force field/database.
- Critical constants for the Peng-Robinson fugacity: `CO2.def`'s header,
  Tc = 304.1282 K, Pc = 7,377,300 Pa, acentric factor = 0.22394.

## Why simulation-vs-simulation, not vs. experiment

Same reasoning as the methane/IRMOF-1 reference (see
`tests/validation/data/irmof1/PROVENANCE.md`): comparing against another
GCMC simulation's output isolates simulation-level discrepancies (sampling,
force field, Ewald convergence) from the additional, harder-to-attribute
gap a real experimental isotherm would introduce (crystal defects,
activation state, force-field transferability).
