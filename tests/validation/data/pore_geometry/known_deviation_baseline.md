# Known-deviation baseline: pore geometry vs. real Zeo++ 0.4.7

Per CLAUDE.md section 4 ("Handling a validated component that disagrees
with its reference"). This file documents where the baseline numbers in
`tests/known_deviation/test_pore_geometry_known_deviation.cc` came from.
See `PROVENANCE.md` in this directory for the full investigation of *why*
each gap exists (not just that it exists) -- in particular the four
independent checks performed on PTY (triclinic lattice vs. Zeo++'s own box
vectors, minimum image vs. 200,000-trial brute force, Voro++ grid
resolution, exact total-tessellation-volume match) before its PLD/ASA/AV
gap was accepted as a genuine, explained methodological difference rather
than debugged further as a suspected bug.

## Provenance

Recorded from a real run of `aleator_known_deviation_tests "[pore]"` (and,
for the tight PLD/ASA/AV numbers referenced for comparison,
`aleator_validation_tests "[pore]"`) under the `dev` preset on 2026-08-06.
LCD is a deterministic, non-Monte-Carlo graph computation, so it needed no
repeat-run check. PTY's ASA/AV baseline uses `PoreAnalysisOptions`'
defaults (`asaSamplesPerAtom=2000`, `volumeSamplesTotal=50000`,
`sampleSeed=0`) -- deterministic per CLAUDE.md invariant #5, so this
reproduces bit-for-bit given the same code.

## LCD (all four structures) -- power (Laguerre) vs. Apollonius diagram

| Structure | computed LCD (Ang) | Zeo++ LCD (Ang) | relative deviation |
|---|---|---|---|
| LTA | 10.26579 | 10.65071 | 3.61% |
| MFI | 5.92922  | 5.95565  | 0.44% |
| FAU | 10.70181 | 10.83869 | 1.26% |
| PTY | 4.67602  | 4.71404  | 0.81% |

LCD is a single global maximum over Voronoi-network free radii. Voro++
(and Zeo++, built on the same underlying approach) computes the power
(Laguerre) diagram for unequal-radius spheres, partitioning space by
`|p-c|^2 - r^2` -- not the "true" Apollonius (additively-weighted) diagram
(`|p-c| - r`) that "largest empty sphere for polydisperse spheres" means
literally, which has curved cell boundaries and isn't what any fast
library computes. Ruled out as a missed-candidate-point bug: also tracking
the maximum free radius seen while sampling edge interiors (already
computed for PLD) never exceeded the vertex-only maximum for LTA. Gap size
does not correlate simply with packing density or cell size across the
four structures tested, consistent with it being a genuine per-structure
geometric-approximation effect rather than a parameter that could be
tuned away.

## PLD / ASA / AV -- LTA, MFI, FAU are tight; PTY is not

| Structure | PLD computed | PLD Zeo++ | PLD dev. | total area computed | total area Zeo++ | area dev. | total volume computed | total volume Zeo++ | volume dev. |
|---|---|---|---|---|---|---|---|---|---|
| LTA | 3.81368 | 3.81365 | 0.00003 Ang (tight) | 203.554 | 192.417 | 5.79% (tight) | 187.510 | 189.203 | 0.90% (tight) |
| MFI | 4.24685 | 4.24727 | 0.00042 Ang (tight) | 349.414 | 359.163 | 2.71% (tight) | 96.826  | 91.614  | 5.69% (tight) |
| FAU | 6.95049 | 6.95049 | 0.00000 Ang (tight) | 1951.56 | 1954.79 | 0.17% (tight) | 2485.79 | 2483.77 | 0.08% (tight) |
| PTY | 3.85082 | 3.92186 | **1.81%** (known deviation) | 50.473 | 45.075 | **11.98%** (known deviation) | 10.7913 | 9.9205 | **8.78%** (known deviation) |

LTA/MFI/FAU are validated TIGHT in `tests/validation/test_pore_geometry.cc`
(PLD within 1e-3 Ang, ASA/AV total within 10%). PTY's numbers are real,
reproducible, and traced to a specific mechanism (not left as an
unexplained gap): the percolation-critical Voronoi edge that sets PTY's
PLD connects two vertices of radius 1.97537 and 1.97491 Ang, but this
codebase's edge-interior sampling (the same machinery validated on LTA)
finds a real constriction of 1.92541 Ang partway along that edge's
physical path -- below both endpoints. Zeo++'s reported threshold,
1.96093 Ang, is suspicious in a specific way: it exactly matches the free
radius of one of this codebase's own Voronoi vertices (node 59). That is
consistent with Zeo++'s percolation graph not fully accounting for this
specific interior-of-edge constriction -- a real methodological
difference in what each tool's percolation graph considers, not a bug in
either. Because PLD is stricter here, the accessible/pocket classification
downstream (ASA/AV) also shifts, which is why PTY's area/volume gaps are
larger than LTA/MFI/FAU's MC-noise-only gaps.

## Not yet independently confirmed

The vertex-radius coincidence for PTY (Zeo++'s reported PLD/2 exactly
matching this codebase's own vertex 59) is strong circumstantial evidence
for the "vertex vs. edge-interior percolation" explanation, not proof --
Zeo++'s own source was not read in this milestone (it was used as a live,
run reference tool). If a future session inspects Zeo++'s percolation
graph construction directly and finds a different explanation, this file
and PROVENANCE.md should be corrected accordingly.
