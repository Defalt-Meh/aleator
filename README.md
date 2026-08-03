# Aleator

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
| `core/geometry` | General triclinic `Lattice`: fractional↔Cartesian conversion, minimum-image convention, cutoff validation, lattice reduction | Brute-force search over up to 125 periodic images on cubic and pathological (60°/60°/60°) triclinic cells |
| `core/neighbor` | `CellList` + `VerletList` with a displacement-based rebuild trigger | Exact agreement with brute-force O(N²) pair search, cubic and triclinic, including a continuous-motion rebuild-trigger stress test |
| `core/math` | Philox4x32-10 counter-based RNG, one independent reproducible stream per (seed, thread) | Bit-exact match against the official Random123 known-answer-test vectors |
| `io` | CIF reader with full space-group symmetry expansion from the asymmetric unit | Real structures from the IZA zeolite database (LTA — cubic, 48 symmetry operations; PTY — triclinic P-1), plus malformed-input error handling |
| `forcefield/pairwise` | Lennard-Jones energy, forces, and virial — truncated / shifted / linear-force-shifted, with analytic long-range tail corrections and Lorentz-Berthelot or geometric mixing rules | Real NIST Standard Reference Simulation Website (SRSW) reference configurations and energies; analytic forces vs. central finite difference |
| `forcefield/electrostatics` | Standard Ewald summation: real-space, reciprocal-space, self-energy, intramolecular exclusion correction, tinfoil boundary conditions | NaCl rock-salt Madelung constant (1.747564594633…, matched to 1.5×10⁻⁸ relative — required 10⁻⁶); NIST SRSW SPC/E water reference energies, term by term; invariance under the Ewald splitting parameter; forces vs. finite difference |
| `engines/monte_carlo` | Grand-canonical Monte Carlo: insertion, deletion, translation, and rotation moves for rigid (single- or multi-site) adsorbate molecules, with a Peng-Robinson equation of state supplying the fugacity in the chemical-potential term | *Validated*: detailed balance verified as a runtime-checked algebraic identity for every move type; Widom test-particle-insertion Henry coefficient matches this engine's own low-pressure isotherm slope to 0.006% on a synthetic system (required 2%), and separately matches this engine's own real IRMOF-1/methane low-pressure loading to within its tight tolerance. *Validated with known deviation*: the full four-point methane/IRMOF-1 isotherm sits systematically 12–15% below the published pyIAST reference curve — see below. |

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

## What's declared but not yet implemented

These have interfaces defined (so the rest of the codebase can be written against them)
but calling them throws `NotImplemented` rather than doing anything:

- **Molecular dynamics** (`engines/dynamics`): velocity-Verlet integration, thermostats,
  barostats. Will be validated on NVE energy drift and equipartition before being trusted.
- **Porous-material geometry analysis** (`engines/geometry_analysis`): Voronoi-based pore
  limiting diameter, largest cavity diameter, accessible surface area. Will be validated
  against Zeo++'s published values for LTA, MFI, and FAU.
- **Structure file writers** (`io`): PDB and LAMMPS `data` output. Reading (CIF) works;
  writing does not yet.
- **Energy-biased Monte Carlo move variants** and multi-species GCMC mixtures.
- **GCMC with electrostatics (charged adsorbates/frameworks).** `Ewald`'s
  reciprocal-space term is a global sum over every particle's structure
  factor and isn't decomposable into a per-particle contribution, so it
  doesn't implement the single-particle trial-move energy GCMC's
  insertion/deletion/translation/rotation moves need. This is enforced at
  `MonteCarloEngine` construction, not discovered mid-run: any force field
  whose `supportsSingleParticleEnergy()` is false (Ewald's default) is
  rejected immediately with an error naming the class. A charged-system
  GCMC isotherm (e.g. CO₂ in Cu-BTC) is therefore not yet possible with
  this engine.
- **MD force-field selection in the CLI's `md` config schema**: `aleator md run` fully
  validates its config and structure file and wires up the integrator, but the
  integrator itself is unimplemented (see above), so a real run always ends in a clean
  `NotImplemented` rather than doing anything — the same is true of `aleator pore
  analyze`, for the same underlying reason.

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
    geometry_analysis/ # pore geometry (declared only)
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
  search) is Highway-vectorized yet. The "faster than RASPA" performance work in CLAUDE.md
  section 5 hasn't started; the SIMD layer today is scaffolding with one real user, not a
  hot-path accelerator.
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

The CLI's design goal is to not repeat RASPA's and LAMMPS's biggest usability problems:
run configuration is plain TOML (never a bespoke scripting language), every config is
fully validated — every required key checked, every value range-checked — before
anything runs, and a malformed config fails immediately with the exact key and line
number, rather than partway through a run that might otherwise take hours.

```bash
aleator --version
aleator gcmc run <config.toml> [--dry-run] [--json]     # grand-canonical Monte Carlo
aleator pore analyze <config.toml> [--dry-run] [--json] # pore geometry (not implemented yet)
aleator md run <config.toml> [--dry-run] [--json]       # molecular dynamics (not implemented yet)
aleator validate <config.toml>                          # validate a config, then exit
aleator bench [--json]                                  # a quick built-in timing check
```

`--dry-run` runs every check a real run would (parsing the config, opening and parsing
the structure file, checking that force-field parameters cover every element actually
present) and prints the fully-resolved configuration, without touching the physics.
`--json` switches the final result to one line of machine-readable JSON on stdout;
progress and diagnostic messages always go to stderr, so piping `--json` output is safe.
Exit codes: `0` success, `1` a configuration or usage error, `2` the requested physics is
genuinely not implemented yet in this build (surfaced as a clean, expected error rather
than a crash or a silent no-op).

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

`examples/` contains real, runnable configs, including a genuine (if deliberately short)
methane-in-IRMOF-1 GCMC run using the same real structure and force field as the
validation suite above — run `aleator gcmc run examples/gcmc_ch4_irmof1.toml` to see it
end to end.

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
  fully explained yet (currently: the CH4/IRMOF-1 isotherm, see above). These never
  turn green by loosening a tolerance around the disagreement; they compare against
  a checked-in baseline and fail only if the gap has grown. CI reports this tier
  separately from `validation` so a documented, tracked disagreement can never mask
  a real regression elsewhere, and a real regression here can never hide inside a
  suite that's "supposed to" have some slack.

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
