add_executable(aleator_known_deviation_tests
    test_gcmc_ch4_irmof1_known_deviation.cc
    test_gcmc_co2_irmof1_known_deviation.cc
    test_pore_geometry_known_deviation.cc
)
aleator_configure_target(aleator_known_deviation_tests)

# Reuses the IRMOF-1 structure/reference data checked in under
# tests/validation/data/ rather than duplicating it.
target_compile_definitions(aleator_known_deviation_tests
    PRIVATE
        ALEATOR_KNOWN_DEVIATION_DATA_DIR="${CMAKE_CURRENT_SOURCE_DIR}/../validation"
)

target_link_libraries(aleator_known_deviation_tests
    PRIVATE
        Catch2::Catch2WithMain
        aleator::core::common
        aleator::core::geometry
        aleator::core::neighbor
        aleator::core::memory
        aleator::core::math
        aleator::engines::monte_carlo
        aleator::engines::geometry_analysis
        aleator::forcefield::electrostatics
        aleator::forcefield::pairwise
        aleator::io
)

include(Catch)
# A separate CTest label from "validation" -- CLAUDE.md section 4: an
# informational check against a checked-in baseline must never be able to
# turn the tight validation suite green by being loose, and CI reports the
# two separately (`ctest -L validation` vs. `ctest -L known-deviation`).
# Deliberately a distinct target/binary rather than a Catch2-tag-derived
# CTest label (ADD_TAGS_AS_LABELS) so this works identically on both the
# vcpkg-pinned Catch2 (3.15.2) and the older Catch2 CMake integration
# shipped by Ubuntu 24.04's system package (3.4.0), which predates
# ADD_TAGS_AS_LABELS (added in Catch2 3.8.0) -- see the "system"
# (non-vcpkg) CMake preset.
catch_discover_tests(aleator_known_deviation_tests
    PROPERTIES LABELS "known-deviation"
    DISCOVERY_MODE PRE_TEST
)
