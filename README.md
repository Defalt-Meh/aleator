# Aleator

A production-grade molecular simulation engine in C++23: Monte Carlo (GCMC/NVT/NPT),
molecular dynamics, and porous-material geometry analysis, sharing one high-performance
periodic-system core. Distributed as a native CLI and as a pip-installable Python package.

See [CLAUDE.md](CLAUDE.md) for the full design contract this project is held to
(unit conventions, validation anchors, portability rules, working agreement).

## Status: scaffolding only — no physics is implemented

This repository currently contains build infrastructure and declared interfaces, not a
working simulation engine. Concretely:

- **Implemented and tested:** SoA particle storage, a triclinic `Lattice` type (matrix
  storage + volume only), the neighbor-list/force-field/counter-based-RNG interfaces,
  an aligned allocator + bump arena, a TOML run-config loader, and one reference Highway
  runtime-SIMD-dispatch kernel (integer vector sum) used as the template every future hot
  kernel will follow.
- **Declared but not implemented — every entry point throws `aleator::NotImplemented`
  rather than a fabricated result:** Lennard-Jones and Ewald energy/forces, the
  Philox counter-based RNG, MD integration, Monte Carlo moves, pore-geometry analysis,
  and CIF/PDB/LAMMPS structure I/O.
- No physics validation test in `tests/validation/` passes a real reference value yet
  (the Madelung-constant, NIST SRSW, and NVE-drift anchors in CLAUDE.md #4 are not wired
  up). The one test currently under that label only asserts that `Ewald::computeEnergy`
  throws instead of returning a plausible-looking number.

Do not use this for scientific work yet.

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

`CMakePresets.json` defines four presets: `dev` (Debug, warnings-as-errors), `release`,
`asan` (AddressSanitizer + UndefinedBehaviorSanitizer, GCC/Clang only), and `bench`
(Release with `-march=native` scoped to the benchmark executable only — local
performance profiling, never shipped; see CLAUDE.md #2.3).

## Building (Python)

```bash
CMAKE_ARGS="-DCMAKE_TOOLCHAIN_FILE=$VCPKG_ROOT/scripts/buildsystems/vcpkg.cmake" \
  pip install -e .
python -c "import aleator; print(aleator.__version__)"
```

The Python extension is built with [nanobind](https://github.com/wjakob/nanobind) via
[scikit-build-core](https://github.com/scikit-build/scikit-build-core). It currently
exposes only `aleator.__version__` and `aleator.vector_sum` (the same reference Highway
kernel mentioned above) — no simulation API yet.

## License

BSD-3-Clause. See [LICENSE](LICENSE).
