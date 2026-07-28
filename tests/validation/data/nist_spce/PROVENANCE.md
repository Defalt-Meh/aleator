# Source

Fetched 2026-07-28 from the NIST Standard Reference Simulation Website (SRSW):

- Page: https://www.nist.gov/mml/csd/chemical-informatics-group/spce-water-reference-calculations-10a-cutoff
- Archive: https://www.nist.gov/document/spcesampleconfigurations-targz

`spce_sample_config_periodic{1,2,3,4}.txt` and `metadata.README` are the
unmodified contents of that archive. The reference energy breakdown
(E_real, E_fourier, E_self, E_intra, E_disp, E_LRC, E_total, all per k_B in
K) these configurations are checked against (in
`tests/validation/test_ewald_spce_water.cc`) is transcribed from the
reference-calculations table on the page above — published, not invented,
per CLAUDE.md invariant #1.

Model parameters (also from that page, SPC/E rigid 3-site water):
q = 0.42380 e (O = -2q, H = +q), sigma_O = 0.316555789 nm,
epsilon_O/k_B = 78.19743111 K, O-H bond = 1 Angstrom, H-O-H angle =
109.47 degrees. LJ site on oxygen only.

Ewald parameters used by NIST for this reference table: alpha =
5.6 / min(Lx, Ly, Lz), kmax = 5 with spherical reciprocal-vector
truncation n0^2+n1^2+n2^2 < kmax^2+2, real-space cutoff = 10.0 Angstrom,
tinfoil (conducting) boundary conditions. NIST's stated physical constants
for this table are CODATA 2010 (k_B = 1.3806488e-23 J/K, e =
1.602176565e-19 C, eps0 = 8.854187817e-12 F/m), a ~1.5e-6 relative
difference from the CODATA 2018/SI-2019 constants
(`Ewald::kCoulombConstant`) used elsewhere in this codebase — accounted
for in the test's tolerance rather than re-deriving a second constant.

NIST's own caveat, transcribed verbatim: "we have reported six significant
figures... the reported energies may differ in the final digits from a
user's calculation" — the test's tolerance is calibrated to that stated
precision limit, not tightened to force a pass (CLAUDE.md invariant #1).
