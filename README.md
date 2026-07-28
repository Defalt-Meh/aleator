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
| `engines/monte_carlo` | Grand-canonical Monte Carlo: insertion, deletion, translation, and rotation moves for rigid (single- or multi-site) adsorbate molecules, with a Peng-Robinson equation of state supplying the fugacity in the chemical-potential term | Detailed balance verified as a runtime-checked algebraic identity for every move type; Widom test-particle-insertion Henry coefficient matches this engine's own low-pressure isotherm slope to 0.006% (required 2%); a full simulated methane adsorption isotherm in IRMOF-1 (real crystal structure, real UFF/TraPPE-derived force field parameters) compared point-by-point against an independently published GCMC-simulated reference isotherm |

The IRMOF-1/methane comparison is the most demanding validation in the codebase and is
reported honestly rather than rounded up: the computed isotherm sits systematically
12–15% below the published curve across the full tested pressure range, a gap too large
to be Monte Carlo sampling noise but consistent with a difference in long-range
dispersion (tail-correction) convention between the two simulations rather than a defect
in the sampling algorithm itself — the sampling machinery is independently confirmed
correct by the exact detailed-balance and Henry-coefficient checks above. See the
comments in `tests/validation/test_gcmc_ch4_irmof1_isotherm.cc` for the full analysis
and the data provenance in `tests/validation/data/irmof1/PROVENANCE.md`.

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
- **CLI subcommands beyond wiring**: the `aleator run` entry point loads a config and
  exercises the engine plumbing end to end, but doesn't yet drive a real GCMC or MD run
  from the command line.
- **Energy-biased Monte Carlo move variants** and multi-species GCMC mixtures.

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
cmake --preset asan && cmake --build --preset asan && ctest --preset asan
cmake --preset bench && cmake --build --preset bench && ./scripts/run_benchmarks.sh
```

`CMakePresets.json` defines four presets:

- `dev` — Debug build, warnings-as-errors. Day-to-day development.
- `release` — optimized build for normal use.
- `asan` — AddressSanitizer + UndefinedBehaviorSanitizer (GCC/Clang only).
- `bench` — Release with `-march=native`, scoped strictly to the local benchmark
  executable. Never used for a distributed binary or wheel; every shipped build uses
  runtime CPU dispatch (via [Highway](https://github.com/google/highway)) instead of a
  compile-time ISA assumption, so the same binary runs correctly — and takes the fastest
  available SIMD path — across SSE4, AVX2, AVX-512, and NEON targets.

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

```bash
aleator --version
aleator run <config.toml>
```

Run configuration is plain TOML rather than a bespoke scripting language:

```toml
[run]
name = "example"
output_directory = "out"
rng_seed = 42
thread_count = 4
```

## Testing

Tests are organized into three tiers, run via CTest with matching labels:

- **unit** — fast, isolated checks of individual types and functions.
- **integration** — cross-module wiring (e.g. that `core/` + `forcefield/` +
  `engines/` actually link and compose correctly).
- **validation** — the physics correctness suite described in the table above: real
  published reference data, analytically known constants, and runtime-checked
  invariants like detailed balance, not just "does it run without crashing."

```bash
ctest --preset dev -L unit
ctest --preset dev -L integration
ctest --preset dev -L validation
```

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
