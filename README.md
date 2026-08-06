# Aleator
## A periodic-system core for Monte Carlo, molecular dynamics, and porous-material geometry analysis 

A molecular simulation engine in modern C++23 covering three domains that are usually
separate tools — grand-canonical and canonical Monte Carlo, molecular dynamics, and
porous-material pore-geometry analysis — built on one shared, high-performance periodic
system core. Distributed as a native command-line binary and as a pip-installable Python
package.

Aleator targets the same problem space as RASPA, LAMMPS, and Zeo++, with a single design
bet: if the periodic-cell geometry, neighbor-list, and force-field core underneath all
three domains is fast and genuinely correct, each simulation engine on top of it becomes
substantially simpler to build and to trust.

## Why this project exists

Molecular simulation results end up in scientific papers, and a wrong number that looks
plausible is worse than an honest error message. Every physics component in this
codebase is built the same way: a validation test against a known, published, or
analytically-derived reference value is written *before* the implementation, and nothing
is considered working until that test passes for real, against real numbers — not
"close enough," not a loosened tolerance. Where a feature isn't implemented yet, calling
it throws a clear `NotImplemented` error rather than silently returning zero or a
plausible-looking guess.

## What's implemented and validated today

| Component | What it does | Validated against |
|---|---|---|
| `core/geometry` | General triclinic `Lattice`: fractional↔Cartesian conversion, minimum-image convention, cutoff validation, lattice reduction, plus an exact O(1) fast path for orthorhombic cells (see Performance below) | Brute-force search over up to 125 periodic images on cubic and pathological (60°/60°/60°) triclinic cells, including cells barely-but-genuinely skewed off 90° (to confirm the fast path is never taken where it shouldn't be) |
| `core/neighbor` | `CellList` + `VerletList` with a displacement-based rebuild trigger | Exact agreement with brute-force O(N²) pair search, cubic and triclinic, including a continuous-motion rebuild-trigger stress test |
| `core/math` | Philox4x32-10 counter-based RNG, one independent reproducible stream per (seed, thread) | Bit-exact match against the official Random123 known-answer-test vectors |
| `io` | CIF reader with full space-group symmetry expansion from the asymmetric unit | Real structures from the IZA zeolite database (LTA — cubic, 48 symmetry operations; PTY — triclinic P-1), plus malformed-input error handling |
| `forcefield/pairwise` | Lennard-Jones energy, forces, and virial — truncated / shifted / linear-force-shifted, with analytic long-range tail corrections and Lorentz-Berthelot or geometric mixing rules | Real NIST Standard Reference Simulation Website (SRSW) reference configurations and energies; analytic forces vs. central finite difference |
| `forcefield/electrostatics` | Standard Ewald summation: real-space, reciprocal-space, self-energy, intramolecular exclusion correction, tinfoil boundary conditions | NaCl rock-salt Madelung constant (1.747564594633…, matched to 1.5×10⁻⁸ relative — required 10⁻⁶); NIST SRSW SPC/E water reference energies, term by term; invariance under the Ewald splitting parameter; forces vs. finite difference |
| `engines/monte_carlo` | Grand-canonical Monte Carlo: insertion, deletion, translation, and rotation moves for rigid (single- or multi-site) adsorbate molecules, with a Peng-Robinson equation of state supplying the fugacity in the chemical-potential term. Supports charged adsorbates/frameworks via an incrementally-maintained Ewald reciprocal-space cache (`EwaldIncrementalState`) — see below. | *Validated*: detailed balance verified as a runtime-checked algebraic identity for every move type; Widom test-particle-insertion Henry coefficient matches this engine's own low-pressure isotherm slope to 0.006% on a synthetic system (required 2%), and separately matches this engine's own real IRMOF-1/methane low-pressure loading to within its tight tolerance. *Validated with known deviation*: the full four-point methane/IRMOF-1 isotherm sits systematically 12–15% below the published pyIAST reference curve — see below. *Validated with known deviation (charged)*: real CO2/IRMOF-1 loading (Ewald electrostatics, real DDEC framework charges) matches a published reference to well under 1σ at low/mid pressure; the high-pressure point deviates by ~2.2σ (~19.5%) — see below. |
| `engines/geometry_analysis` | Porous-material geometry: largest cavity diameter (LCD), pore limiting diameter (PLD), N₂-probe accessible surface area (ASA) and accessible volume (AV), from a periodic radical (power/Laguerre) Voronoi decomposition of the framework via Voro++, probe-radius parameterized, triclinic cells first-class. | *Validated*: PLD and ASA/AV match real Zeo++ 0.4.7 output on real IZA zeolite structures (LTA, MFI, FAU) — PLD to ≤4×10⁻⁴ Å, ASA/AV within 10%; `Lattice::minimumImageDisplacement` cross-checked against a 200,000-trial brute-force periodic-image search on a genuinely triclinic structure (PTY); FAU's real published inaccessible sodalite-cage pockets correctly come out as a nonzero, minority fraction of the void space. *Validated with known deviation*: LCD sits systematically 0.4–3.6% below Zeo++ on every structure tested (power-vs-Apollonius diagram, a real algorithmic approximation shared with Zeo++'s own tessellation library); PTY's PLD/ASA/AV additionally sit 1.8%/12%/8.8% off Zeo++, traced to a specific percolation-critical Voronoi edge, not left unexplained — see below. |

The IRMOF-1/methane comparison against the published pyIAST reference isotherm is the
most demanding validation in the codebase and is reported honestly rather than rounded
up. Per CLAUDE.md section 4, it is split into two separate CTest-visible tests: a tight
test (`tests/validation/test_gcmc_ch4_irmof1_isotherm.cc`) checking this engine's own
internal self-consistency on the real structure (no external reference involved, so a
regression in the sampling machinery has nowhere to hide behind a wide tolerance), and
an informational `known-deviation` test
(`tests/known_deviation/test_gcmc_ch4_irmof1_known_deviation.cc`) that compares the full
curve against pyIAST's published numbers and fails only if the gap widens beyond a
checked-in baseline — never by being loose about the gap itself.

The gap is real, systematic (12–15% below reference at every one of four pressure points
spanning a 100x range), and **unresolved**. A same-session investigation ruled out
several standing hypotheses rather than assuming them: RASPA2's own GenericMOFs
force-field file (re-fetched live) states "shifted" truncation with "no tail
corrections" — exactly what this codebase already does, so a *missing* tail correction
is not an available explanation; the Lennard-Jones parameters and Lorentz-Berthelot
mixing rule were re-verified against that same live-fetched source and match exactly;
inaccessible-pore blocking was ruled out by reading this codebase's own insertion code
(uniform sampling over the full cell, no accessibility filtering, so blocking can only
affect sampling efficiency, not the equilibrium loading); and a real sensitivity check
(cutoff 12.0 vs. 12.8 Å at constant cell size) moved the computed loading by only 0.25
standard errors, ruling out the within-cell portion of the cutoff-margin hypothesis. A
larger cutoff enabled by a bigger simulation cell, and pyIAST's own (publicly
undocumented) simulation parameters, remain open and unverified. Full writeup:
`tests/validation/data/irmof1/known_deviation_baseline.md`.

### Charged adsorbates: incremental Ewald for GCMC

The reciprocal-space Ewald sum is a global function of every charge's structure factor
S(k) — it does not decompose into a per-particle contribution the way the real-space
term does, so a GCMC trial move naively recomputing it from scratch would cost O(N) per
trial (defeating the point of trial moves) or worse. `forcefield/electrostatics/
ewald_incremental_state.hpp` solves this with a dedicated, engine-owned cache: S(k), the
self energy, and the intramolecular exclusion correction are maintained incrementally —
`propose*()` computes a trial energy delta without mutating anything, `commit*()` folds
it in only on acceptance, so a rejected move needs no rollback logic at all (there is
nothing to undo). This is deliberately **not** exposed through `ForceField::
computeParticleEnergy()`/`supportsSingleParticleEnergy()`: that interface's contract is
"the complete per-particle energy," which the reciprocal term structurally cannot supply
that way, so `Ewald` honestly keeps advertising `false` there and `MonteCarloEngine`
gained a separate, explicitly-typed optional `Ewald` component instead. Real-space
Coulomb, in contrast, genuinely is pairwise/local, so it stays a direct O(N) scan
(`Ewald::realSpaceParticleEnergy`), the same cost structure as the dispersion force
field.

Validated: a real charged system (the same real IRMOF-1 framework, now carrying its
real DDEC framework charges, plus rigid CO2 molecules) run through `MonteCarloEngine`
for 10⁵ real GCMC steps agrees with a full from-scratch recomputation of the cached
energy to machine precision throughout (drift 0.00×10⁰–1.79×10⁻¹⁶ relative, gated at
1×10⁻¹⁰); a rejected trial leaves the cache provably bit-identical; and inserting then
deleting the same molecule produces exactly opposite energy deltas (the energy-level
identity the detailed-balance ratio proof depends on). Long-run floating-point
accumulation in the incremental cache is bounded by periodic resync against a fresh
recomputation — every 100 accepted moves in production, an interval set from a real
measurement (500 was enough at a small synthetic system's scale but not at the real
~850-charge IRMOF-1/CO2 system's scale, where it left ~9×10⁻¹⁰ relative drift).

The real published-reference test — CO2 in IRMOF-1 at 298 K, since HKUST-1/Cu-BTC (this
milestone's other candidate system) was checked first and confirmed absent from the only
available real GCMC-simulated reference dataset — reproduces the same "validated with
known deviation" honesty as the methane case: computed loading at 10⁴ and 10⁵ Pa matches
the published reference (the CRAFTED database's real RASPA GCMC simulation of this exact
system) to well under 1 combined standard error (~0.09σ and ~0.19σ); at 10⁶ Pa
(~127 adsorbed molecules at equilibrium) the computed loading sits ~19.5% (~2.2σ) below
the reference. The gap is tracked, not hidden: the charged-GCMC energy machinery itself
is independently verified correct (the drift gate above), and the leading suspected
cause — this test's equilibration length being calibrated for low occupancy and likely
insufficient at ~127 molecules — is disclosed as unconfirmed rather than assumed. Full
writeup: `tests/validation/data/co2_irmof1/known_deviation_baseline.md`.

### Performance: the CH4/IRMOF-1 isotherm, 32 s → 3.97 s

The four-point CH4/IRMOF-1 known-deviation isotherm
(`tests/known_deviation/test_gcmc_ch4_irmof1_known_deviation.cc`, 480,000 total GCMC
steps) is this codebase's one already-validated headline workload, and "faster than
RASPA" is the pitch (CLAUDE.md section 5), so it is the one benchmarked here.
Before this milestone it took 32.0 s (Release) / 58 s (Debug, `dev` preset). Profiling
first (`sample`, macOS; a real 25 s stack-sample window, not a guess) showed **91.2%
of all time inside `Lattice::minimumImageDisplacement`** — called from a direct
O(frameworkCount) scan on every trial move, `LennardJones::pairEnergy` itself was only
6.8%. Two things were tried; only the second one actually moved the needle for this
system:

1. **Guest-guest interactions wired through `core::CellList`** (translation/insertion/
   deletion now query a `CellList` built over the mobile guest-molecule tail of
   `particles_` instead of scanning every particle). Proven an *exact* refactor, not an
   approximation: `LennardJones::computeParticleEnergyOverCandidates` sums over a
   candidate list instead of `0..N`, which is bit-identical to the old full scan as long
   as (a) every candidate beyond `cutoff()` contributes exactly `0.0` (true by
   construction — IEEE754 `x + 0.0 == x` exactly, so dropping them changes nothing) and
   (b) every candidate that *does* contribute is visited in the same relative order the
   old ascending scan visited it (floating-point addition is not associative, so
   `core::CellList`'s bin-traversal order is explicitly sorted back to index order
   before use). Verified directly (`tests/validation/test_gcmc_cell_list_bit_identical.cc`,
   `REQUIRE(a == b)`, not a tolerance, across cubic and triclinic cells and guest
   occupancy up to CO2/IRMOF-1's real ~127-molecule scale) and end-to-end (the whole
   known-deviation test reproduces bit-identical loading values before and after, same
   seed). On *this* test, though, it made **no measurable difference** — guest occupancy
   here tops out around 23 molecules, too few for an O(N)-vs-cell-list distinction to
   show up against a 424-atom framework scan. It matters for the higher-occupancy
   charged CO2/IRMOF-1 system and is real, validated infrastructure either way (CLAUDE.md
   section 3's "cash the bet" on `core/neighbor` actually being used by an engine).

2. **A rigid-framework guest-host energy grid was built, validated, and *not* adopted**
   (`engines/monte_carlo/framework_energy_grid.hpp`, `FrameworkEnergyGrid`): one
   trilinearly-interpolated energy table per guest LJ species, replacing the
   O(frameworkCount) direct scan with an O(1) lookup — the literal "tabulate the guest-
   host field once" idea. It works (a real accuracy-vs-spacing-vs-build-time sweep is in
   `tests/validation/test_framework_energy_grid.cc`, with a raw-node cap
   (`kEnergyCapKelvin = 1e5`) needed because a grid node landing near a framework atom's
   core computes a real 1e17+ K raw LJ repulsion that a naive interpolation smears into
   meaningless output — physically harmless to cap, since `exp(-1e5 / 298 K)` is already
   indistinguishable from zero acceptance probability). But **measured, not assumed**: at
   this exact 25.832 Å single-unit-cell system (12.0 Å cutoff, close to `L_perp/2`), a
   spacing coarse enough to build quickly (0.5 Å, ~11 s) was *not* accurate enough — it
   more than doubled this test's own tight self-consistency deviation past its 12% gate
   (22.6%, vs. 3.4% without it) — and a spacing accurate enough to pass (0.2 Å) cost
   **~165 s to build once**, over 5× this entire 4-point isotherm's *pre-optimization*
   runtime. Root cause: a single unit cell with `cutoff ≈ L_perp/2` means nearly the
   whole framework is "in range" of any point in the cell, so neither a cell list nor a
   precomputed grid can avoid touching most of it — this is the same underlying
   constraint the methane/IRMOF-1 known-deviation writeup and the CO2/IRMOF-1 supercell
   already ran into (see `CLAUDE.md` section 0, defect 1). `FrameworkEnergyGrid` is real,
   tested infrastructure (and `engines/monte_carlo/monte_carlo_engine.hpp`'s GCMC
   constructor accepts one), scoped for systems where `cutoff << L` or where far more
   total steps amortize the one-time build — not wired into this test's default path.

3. **The actual fix, found from following the profile rather than the two techniques
   above: `Lattice::minimumImageDisplacement` had a hidden orthorhombic fast path
   available and wasn't using it.** IRMOF-1's cell is cubic (90°/90°/90°), but the
   general-purpose triclinic algorithm (a 125-candidate search around a Gauss-reduced
   basis — necessary for a genuinely skewed cell, see CLAUDE.md's Ewald/triclinic
   warnings) was running unconditionally, on every call, for every cell shape. An
   orthorhombic cell's minimum image is *exact* via independent per-axis rounding — no
   search needed, because the axes are mutually orthogonal; this does not generalize to
   triclinic cells, which is exactly why the general path still exists and is still used
   for non-orthorhombic input. Detecting "orthorhombic" needs a tolerance, not an exact
   zero check: a real CIF-derived 90°/90°/90° cell's off-diagonal matrix entries are
   never exactly `0.0` (`cos(π/2) ≈ 6.1×10⁻¹⁷` in double precision, not `0`), so the
   check is relative-tolerance-based (`1e-9`, chosen to comfortably absorb that
   trig-rounding noise while staying ~1e5x below the off-diagonal magnitude a cell even
   0.01° off orthorhombic would have) — validated explicitly on a genuinely-if-barely
   skewed (89.5°) cell to confirm it still takes, and is still correct on, the general
   path (`tests/validation/test_lattice_minimum_image.cc`).

**Result**: the fast path alone took the 4-point isotherm from 32.0 s to **3.97 s**
(Release; 58 s → 37.1 s Debug) — under the 5 s target — with **bit-identical** computed
loading values before and after (same seed, same output to every printed digit): this is
an exact optimization, not an approximation, so it needed no known-deviation-baseline
re-justification. `benchmarks/bench_gcmc_ch4_irmof1.cc` tracks a synthetic-system
equivalent of this workload as an ongoing CI-checked baseline (`scripts/
run_benchmarks.sh`, CLAUDE.md's >5% regression gate).

### Pore geometry: periodic radical Voronoi decomposition

`engines/geometry_analysis` computes LCD, PLD, ASA, and AV from a periodic Voronoi
network of the framework, weighted by each atom's van der Waals radius (a "radical" or
power-diagram decomposition, since atoms have different sizes) — the same conceptual
approach Zeo++ itself uses, and the standard building block RASPA/Zeo++ users expect.

**Dependency choice, made and justified before writing code**: [Voro++](https://math.lbl.gov/voro++/)
(BSD-3-Clause-LBNL), not a reimplementation. It's mature, supports periodic triclinic
domains with per-particle radii (`container_periodic_poly`) out of the box, and is
available both via vcpkg and as an Ubuntu 24.04 system package (`libvoro++-dev`) —
consistent with this codebase's "no hard vcpkg requirement" portability contract.

**Triclinic support, made real, not assumed**: Voro++ requires its periodic box in a
canonical lower-triangular form (`(bx,0,0),(bxy,by,0),(bxz,byz,bz)`); this codebase's own
CIF-derived lattice matrix already happens to satisfy that form for the real structures
tested, but the general transform (a Gram-Schmidt orthogonalization of the lattice's own
vectors — an isometric re-expression, not an approximation, since it preserves the
lattice's Gram/metric matrix exactly) is implemented and used unconditionally, with every
position round-tripped through fractional coordinates rather than an explicit rotation
matrix. `Lattice::minimumImageDisplacement` itself — the actual PBC primitive invariant 6
requires be tested on a non-orthogonal cell — is cross-checked against a 200,000-trial
brute-force periodic-image search (-3..3 in each direction) on the real triclinic PTY
structure (α=84.6°, β=83.8°, γ=86.7°): zero mismatches, max relative error 1.5×10⁻¹⁴
(floating-point noise).

**A real bug found and fixed during validation, not assumed away**: an early
implementation weighted each Voronoi-network edge by `min(endpoint radii)`. On real LTA,
this gave a PLD of 10.016 Å (≈ the LCD, i.e. "no real constriction anywhere" — obviously
wrong for a zeolite). Traced to a specific oxygen atom at a special symmetry position
(x=0) whose Voronoi cell has a face with two vertices 11.6 Å apart: the true constriction
is in the *interior* of that long edge, not at either endpoint, and both endpoints happen
to sit in wide cage-center regions on either side of it. Fixed by sampling the free
radius along each edge's actual physical path (spacing chosen empirically: 0.2 Å already
matched Zeo++ to ≤3×10⁻³ Å on LTA/MFI/FAU; the shipped 0.02 Å matches to ≤4×10⁻⁴ Å).

**The honest remainder**: LCD sits 0.4–3.6% below real Zeo++ 0.4.7 output on every
structure tested — a genuine, understood gap between Voro++'s power (Laguerre) diagram
and the "true" Apollonius (additively-weighted) diagram that "Voronoi decomposition of
unequal-radius spheres" technically means (curved cell boundaries, not implemented by any
fast library, including the one underneath Zeo++ itself). And PTY (the triclinic
structure) has a larger PLD/ASA/AV gap than LTA/MFI/FAU, investigated specifically rather
than assumed to be more of the same LCD-style noise: four independent checks (lattice vs.
Zeo++'s own printed box vectors, minimum image vs. brute force, Voro++ grid-resolution
insensitivity, exact total-tessellation-volume match) ruled out a bug before the gap was
traced to a specific percolation-critical edge whose interior constriction this
codebase's finer sampling finds and Zeo++'s reported number appears not to — full writeup
in `tests/validation/data/pore_geometry/PROVENANCE.md` and
`known_deviation_baseline.md`.

## What's declared but not yet implemented

These have interfaces defined (so the rest of the codebase can be written against them)
but calling them throws `NotImplemented` rather than doing anything:

- **Molecular dynamics** (`engines/dynamics`): velocity-Verlet integration, thermostats,
  barostats. Will be validated on NVE energy drift and equipartition before being trusted.
- **Structure file writers** (`io`): PDB and LAMMPS `data` output. Reading (CIF) works;
  writing does not yet.
- **Energy-biased Monte Carlo move variants** and multi-species GCMC mixtures.
- **`[md]` config schema exists (`io/config.hpp`) with no CLI subcommand to run it**:
  `engines/dynamics` is itself unimplemented (see above), and CLAUDE.md's CLI milestone
  is explicit that a subcommand which just prints `NotImplemented` when run is worse than
  not having it — so there is no `aleator md run`. `aleator validate` still
  schema-validates `[md]` config files (real, useful groundwork ahead of that engine
  existing) and says plainly that there's no engine for it yet, distinct from a malformed
  config.

## Architecture

```
src/
  core/          # shared substrate: lattice/PBC, neighbor lists, SoA particle storage,
                  # counter-based RNG, aligned memory, runtime-dispatched SIMD kernels
  forcefield/
    pairwise/          # Lennard-Jones
    electrostatics/    # Ewald summation
    parameters/        # force-field parameter file loading
  engines/
    monte_carlo/       # GCMC moves, acceptance criteria, Peng-Robinson EOS
    dynamics/          # MD integrators (declared only)
    geometry_analysis/ # pore geometry: Voro++-based periodic radical Voronoi network,
                        # percolation graph, LCD/PLD/ASA/AV
  io/            # CIF reading, TOML run configuration
  cli/           # command-line entry point
  bindings/      # Python extension module (nanobind)
```

The dependency direction is enforced by the build: `core/` depends on nothing else in
this repository; `forcefield/` depends only on `core/`; `engines/` depend on `core/` and
`forcefield/` but never on each other laterally; `io/` and `cli/` depend downward only.

### Internal unit convention

One convention, enforced at every I/O boundary rather than left implicit: length in
Ångström, energy in Kelvin (i.e. energy divided by the Boltzmann constant), mass in
atomic mass units, time in picoseconds, charge in elementary charge units. File formats
and external APIs that use different units (CIF's Ångström cells, Pascal-based
equations of state, SI physical constants) convert explicitly at the point they enter or
leave this convention — never implicitly, and never by assuming a caller already knows.

### Determinism

Random numbers come from Philox4x32-10, a counter-based generator: a stream is
identified entirely by a `(seed, streamIndex)` pair rather than mutable shared state, so
independent, reproducible streams per thread or per Monte Carlo trial come for free.
Given the same seed, the same input, and the same thread count, a simulation reproduces
bit-identical output — a property that matters for debugging and for reproducing a
published result exactly.

## Building (C++)

Requires a [vcpkg](https://github.com/microsoft/vcpkg) checkout with `VCPKG_ROOT` set,
CMake ≥ 3.25, and Ninja.

```bash
cmake --preset dev && cmake --build --preset dev
ctest --preset dev                    # all tests
ctest --preset dev -L validation      # physics validation only
ctest --preset dev -L known-deviation # documented, tracked disagreements with published references
cmake --preset asan && cmake --build --preset asan && ctest --preset asan
cmake --preset bench && cmake --build --preset bench && ./scripts/run_benchmarks.sh
```

`CMakePresets.json` defines five presets:

- `dev` — Debug build, warnings-as-errors. Day-to-day development.
- `release` — optimized build for normal use.
- `asan` — AddressSanitizer + UndefinedBehaviorSanitizer (GCC/Clang only).
- `bench` — Release with `-march=native`, scoped strictly to the local benchmark
  executable. Never used for a distributed binary or wheel; every shipped build uses
  runtime CPU dispatch (via [Highway](https://github.com/google/highway)) instead of a
  compile-time ISA assumption, so the same binary runs correctly — and takes the fastest
  available SIMD path — across SSE4, AVX2, AVX-512, and NEON targets. As of this writing
  that runtime dispatch is wired up for exactly one kernel (`vectorSum`, an internal
  integer reduction used in tests) — no physics kernel (Lennard-Jones, Ewald, neighbor
  search) is Highway-vectorized yet; the SIMD layer today is scaffolding with one real
  user, not a hot-path accelerator. CLAUDE.md section 5's "faster than RASPA" performance
  work has started (see the Performance section above) but hasn't touched SIMD — its one
  win so far (the CH4/IRMOF-1 isotherm, 32.0 s → 3.97 s) came from algorithmic fixes
  (an exact orthorhombic minimum-image fast path, a guest-guest cell list), not
  vectorization.
- `system` — same as `dev`, but configured against system-installed dependencies instead
  of vcpkg (no `VCPKG_ROOT` needed; Python bindings off). vcpkg is the supported,
  CI-covered path, but `find_package()` for every dependency also succeeds against a
  system install (e.g. Ubuntu 24.04: `apt install libhwy-dev catch2 libspdlog-dev
  libbenchmark-dev libtomlplusplus-dev`) — a contributor who can't or doesn't want to
  bootstrap vcpkg can still build and run the full test suite. See the
  `system-packages` CI job for the exact, verified package list.

## Building (Python)

```bash
CMAKE_ARGS="-DCMAKE_TOOLCHAIN_FILE=$VCPKG_ROOT/scripts/buildsystems/vcpkg.cmake" \
  pip install -e .
python -c "import aleator; print(aleator.__version__)"
```

The Python extension is built with [nanobind](https://github.com/wjakob/nanobind) via
[scikit-build-core](https://github.com/scikit-build/scikit-build-core), and packaged as
portable wheels (Linux x86-64/aarch64, macOS x86-64/arm64, Windows x86-64) via
`cibuildwheel`. The bound API surface currently mirrors the C++ scaffolding rather than
the full simulation engine — a fuller, ergonomic Python API is planned once the engine
surface above stabilizes further.

## Command-line interface

The CLI is meant to be a real differentiator, not an afterthought: both RASPA and
LAMMPS have genuinely painful UX, and that's a defensible opening. The design goal is
to not repeat their biggest usability problems: run configuration is plain TOML (never
a bespoke scripting language), every config is fully validated — every required key
checked, every value range-checked — before anything runs, a malformed config fails
immediately with the exact key and line number rather than partway through a run that
might otherwise take hours, and **only subcommands backed by a working engine exist at
all**:

```bash
aleator --version
aleator gcmc run <config.toml> [--dry-run] [--json]     # grand-canonical Monte Carlo
aleator pore analyze <config.toml> [--dry-run] [--json] # LCD, PLD, ASA, AV
aleator validate <config.toml>                          # validate a config, then exit
aleator bench [--json]                                  # a quick built-in timing check
```

GCMC and pore-geometry analysis are the implemented engines right now (see the status
table above), so those are the runnable subcommands — there is no `aleator md run`
(`engines/dynamics` is still `NotImplemented`). A subcommand that just prints "not
implemented" when run is worse than not having it: it advertises capability the tool
doesn't actually have. The `[md]` config *schema* still exists (real groundwork for when
that engine lands), and `aleator validate` still recognizes it — it will schema-validate
the file and then say plainly that there's no engine for it yet (exit code `2`, distinct
from a malformed config's `1`), never silently call it "valid" in a way that could be
mistaken for "runnable."

`--dry-run` runs every check a real run would (parsing the config, opening and parsing
the structure file, checking that force-field parameters cover every element actually
present — for `pore analyze`, that every element in the CIF has a known van der Waals
radius) and prints the fully-resolved configuration — every default filled in, not just
the keys the config set explicitly — without touching the physics or writing anything.
`--json` switches that to one line of machine-readable JSON on stdout; progress and
diagnostic messages always go to stderr, so piping `--json` output is safe. A real
(non-dry-run) `gcmc run` or `pore analyze` writes that same fully-resolved configuration,
including the RNG seed, to `<run.output_directory>/resolved_config.json` (resolved
relative to the config file's own location, not whatever directory `aleator` happened to
be invoked from) — so a published result can be reproduced from its artifacts alone. Exit
codes: `0` success, `1` a configuration or usage error, `2` `validate` was pointed at a
schema-valid config for an engine this build doesn't implement yet.

A minimal GCMC config, with an explicit per-element Lennard-Jones parameter table rather
than a hidden built-in force-field database:

```toml
[run]
name = "example"
rng_seed = 42

[gcmc]
framework_cif = "IRMOF-1.cif"
temperature_kelvin = 298.0
pressure_bar = 1.0

[gcmc.adsorbate]
name = "CH4"
epsilon_kelvin = 158.5
sigma_angstrom = 3.72
mass_amu = 16.04246

[[gcmc.framework_lj]]
element = "Zn"
epsilon_kelvin = 62.3992
sigma_angstrom = 2.46155
# ...one [[gcmc.framework_lj]] entry per element actually present in the CIF
```

Every physical quantity in the schema carries its unit in the key name
(`temperature_kelvin`, `pressure_bar`, `cutoff_angstrom`, `epsilon_kelvin`,
`sigma_angstrom`, `mass_amu`) — never a bare `temperature` or `pressure`. The optional
`gcmc.energy_grid_spacing_angstrom` key (CLAUDE.md section 5's `FrameworkEnergyGrid`) is
unset by default, meaning "direct O(N) guest-host scan," not an implicit spacing —
setting it is a real accuracy/speed tradeoff, not free, so it's opt-in (see that
section's README writeup for the measured cost).

`examples/` contains real, runnable configs. `examples/gcmc_ch4_irmof1.toml` is not a
shortened smoke test: it uses the exact same seed, pressure, and step counts as
`tests/validation/test_gcmc_ch4_irmof1_isotherm.cc`'s real 0.1 bar point, and the CLI
samples its running mean the same way that test does (every step, not once per progress
update) — so `aleator gcmc run examples/gcmc_ch4_irmof1.toml` reproduces that test's `<N>`
and loading numbers exactly, not approximately (same deterministic RNG stream, CLAUDE.md
invariant #5), in about a second thanks to CLAUDE.md section 5's performance milestone.
`tests/integration/test_cli_end_to_end.cc` checks this cross-consistency directly.
`examples/pore_irmof1.toml` runs `pore analyze` on the same real IRMOF-1 structure (424
framework atoms after symmetry expansion) — a genuine multi-minute computation at this
size, unlike the GCMC example, so use `--dry-run` there for a quick config check.

## Testing

Tests are organized into four tiers, run via CTest with matching labels:

- **unit** — fast, isolated checks of individual types and functions.
- **integration** — cross-module wiring (e.g. that `core/` + `forcefield/` +
  `engines/` actually link and compose correctly).
- **validation** — the physics correctness suite described in the table above: real
  published reference data, analytically known constants, and runtime-checked
  invariants like detailed balance, not just "does it run without crashing." Every
  entry here is tight: a real regression of ~10% or more will fail it.
- **known-deviation** — informational tests for a validated component that still
  disagrees with an external published reference in a way this codebase hasn't
  fully explained yet (currently: the CH4/IRMOF-1 and CO2/IRMOF-1 isotherms, and
  pore-geometry LCD/PTY, see above). These never turn green by loosening a tolerance
  around the disagreement; they compare against a checked-in baseline and fail only
  if the gap has grown. CI reports this tier separately from `validation` so a
  documented, tracked disagreement can never mask a real regression elsewhere, and a
  real regression here can never hide inside a suite that's "supposed to" have some
  slack.

```bash
ctest --preset dev -L unit
ctest --preset dev -L integration
ctest --preset dev -L validation
ctest --preset dev -L known-deviation
```

`scripts/check_test_discovery.py` (run in CI right after the build, before any test
tier) independently counts every `TEST_CASE` under `tests/` and compares it against
how many CTest actually discovered, failing the build on any mismatch — the guard
against a test silently becoming invisible to the runner (CLAUDE.md invariant #8).

Reference data used by the validation suite (NIST SRSW configurations, IZA zeolite
structures, the IRMOF-1 crystal structure and force field, published isotherm data) is
committed under `tests/validation/data/`, each with a `PROVENANCE.md` recording exactly
where it came from and when it was fetched — nothing in the validation suite is invented
or hand-tuned to pass.

## Portability

Targets Linux (x86-64, aarch64), macOS (x86-64, Apple Silicon), and Windows (x86-64), on
GCC ≥ 13, Clang ≥ 16, and MSVC ≥ 19.38. No POSIX-only APIs, no MPI in the distributed
package (threading only; cluster builds can opt into MPI separately), and fixed-width
integer types throughout rather than assumptions about platform-dependent widths.

## License

BSD-3-Clause. See [LICENSE](LICENSE).

## Citation

See [CITATION.cff](CITATION.cff).
