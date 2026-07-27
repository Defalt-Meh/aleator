# Source

Fetched 2026-07-27 from the NIST Standard Reference Simulation Website (SRSW):

- Page: https://www.nist.gov/mml/csd/chemical-informatics-group/lennard-jones-fluid-reference-calculations-cuboid-cell
- Archive: https://www.nist.gov/document/ljsampleconfigurations-targz

`lj_sample_config_periodic{1,2,3,4}.txt` and `metadata.README` are the
unmodified contents of that archive. The reference energy/virial/tail-
correction values these configurations are checked against (in
`tests/validation/test_lennard_jones_nist.cc`) are transcribed from the
"Reference Calculations for the Internal Energy and Pair Virial" table on
the page above — published, not invented, per CLAUDE.md invariant #1.
