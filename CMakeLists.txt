cmake_minimum_required(VERSION 3.25)

project(aleator
    VERSION 0.1.0
    DESCRIPTION "Monte Carlo and porous-material geometry analysis on a shared periodic-system core"
    HOMEPAGE_URL "https://github.com/Defalt-Meh/aleator"
    LANGUAGES CXX
)

if(PROJECT_SOURCE_DIR STREQUAL PROJECT_BINARY_DIR)
    message(FATAL_ERROR "In-source builds are not allowed. Use a preset, e.g. `cmake --preset dev`.")
endif()

set(CMAKE_CXX_STANDARD 23)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake")

# ---------------------------------------------------------------------------
# Options
# ---------------------------------------------------------------------------
option(ALEATOR_ENABLE_MPI "Enable MPI-based distributed parallelism (cluster builds only, never in the wheel)" OFF)
option(ALEATOR_BUILD_PYTHON "Build the nanobind Python extension module" ON)
option(ALEATOR_BUILD_BENCHMARKS "Build the Google Benchmark suite" ON)
option(ALEATOR_BUILD_TESTS "Build the Catch2 test suite" ON)
option(ALEATOR_BENCH_NATIVE
    "Compile ONLY the benchmarks/ executable with -march=native. Local-only: \
never set for a distributable/wheel/release build. See CLAUDE.md #2.3 and #5."
    OFF
)

if(CMAKE_CXX_FLAGS MATCHES "-march=native")
    message(FATAL_ERROR
        "-march=native was found in the global CMAKE_CXX_FLAGS. It is forbidden "
        "outside the benchmarks/ target (see CLAUDE.md #2.3). Use "
        "-DALEATOR_BENCH_NATIVE=ON instead, which scopes the flag to the "
        "benchmark executable only. Hot kernels must use Highway's runtime "
        "CPU dispatch, never a compile-time ISA assumption."
    )
endif()

if(ALEATOR_BENCH_NATIVE AND NOT ALEATOR_BUILD_BENCHMARKS)
    message(FATAL_ERROR "ALEATOR_BENCH_NATIVE requires ALEATOR_BUILD_BENCHMARKS=ON.")
endif()

include(Warnings)
include(Sanitizers)
include(AleatorTarget)

# ---------------------------------------------------------------------------
# Generated version header (aleator/version.hpp), single source of truth
# taken from project(... VERSION ...) above.
# ---------------------------------------------------------------------------
configure_file(
    ${PROJECT_SOURCE_DIR}/cmake/version.hpp.in
    ${PROJECT_BINARY_DIR}/generated/aleator/version.hpp
    @ONLY
)
add_library(aleator_version INTERFACE)
add_library(aleator::version ALIAS aleator_version)
target_include_directories(aleator_version INTERFACE
    $<BUILD_INTERFACE:${PROJECT_BINARY_DIR}/generated>
    $<INSTALL_INTERFACE:include>
)

# ---------------------------------------------------------------------------
# Dependencies (vcpkg manifest mode: see vcpkg.json)
# ---------------------------------------------------------------------------
find_package(hwy CONFIG REQUIRED)
find_package(tomlplusplus CONFIG REQUIRED)
find_package(spdlog CONFIG REQUIRED)
# engines/geometry_analysis milestone: Voro++ (github.com/chr1shr/voro), the
# standard building block for periodic Voronoi decomposition -- a mature,
# BSD-3-licensed library rather than a reimplementation, since a correct,
# numerically-robust periodic *weighted* (radical/Laguerre, needed because
# framework atoms have different radii) Voronoi decomposition is itself a
# research-grade computational-geometry undertaking. Available both via
# vcpkg (port "voro") and system packages (Ubuntu 24.04: libvoro++-dev, in
# universe) -- verified against the real package index, not assumed, per
# CLAUDE.md's portability contract.
find_package(VORO CONFIG REQUIRED)

if(ALEATOR_ENABLE_MPI)
    find_package(MPI REQUIRED COMPONENTS CXX)
endif()

# ---------------------------------------------------------------------------
# Library sources
#
# Dependency direction (CLAUDE.md #3): core/ depends on nothing above it;
# forcefield/ depends only on core/; engines/ depend on core/ + forcefield/;
# io/ and cli/ depend downward only. Each subdirectory below enforces this
# through the PUBLIC/PRIVATE visibility of its target_link_libraries() call
# — do not add a link that violates it without updating that rule in
# CLAUDE.md first.
# ---------------------------------------------------------------------------
add_subdirectory(src/core)
add_subdirectory(src/forcefield)
add_subdirectory(src/engines)
add_subdirectory(src/io)
add_subdirectory(src/cli)

if(ALEATOR_BUILD_PYTHON)
    add_subdirectory(src/bindings)
endif()

# ---------------------------------------------------------------------------
# Tests / benchmarks
# ---------------------------------------------------------------------------
include(CTest)
if(ALEATOR_BUILD_TESTS AND BUILD_TESTING)
    find_package(Catch2 3 CONFIG REQUIRED)
    add_subdirectory(tests)
endif()

if(ALEATOR_BUILD_BENCHMARKS)
    find_package(benchmark CONFIG REQUIRED)
    add_subdirectory(benchmarks)
endif()

# ---------------------------------------------------------------------------
# Status summary
# ---------------------------------------------------------------------------
message(STATUS "aleator ${PROJECT_VERSION}")
message(STATUS "  C++ standard          : ${CMAKE_CXX_STANDARD}")
message(STATUS "  Build type            : ${CMAKE_BUILD_TYPE}")
message(STATUS "  ALEATOR_ENABLE_MPI    : ${ALEATOR_ENABLE_MPI}")
message(STATUS "  ALEATOR_BUILD_PYTHON  : ${ALEATOR_BUILD_PYTHON}")
message(STATUS "  ALEATOR_BUILD_BENCHMARKS : ${ALEATOR_BUILD_BENCHMARKS}")
message(STATUS "  ALEATOR_BUILD_TESTS   : ${ALEATOR_BUILD_TESTS}")
message(STATUS "  ALEATOR_BENCH_NATIVE  : ${ALEATOR_BENCH_NATIVE}")
message(STATUS "  ALEATOR_ENABLE_ASAN   : ${ALEATOR_ENABLE_ASAN}")
message(STATUS "  ALEATOR_ENABLE_UBSAN  : ${ALEATOR_ENABLE_UBSAN}")
