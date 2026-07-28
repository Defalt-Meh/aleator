add_executable(aleator_validation_tests
    test_cif_real_structures.cc
    test_ewald_alpha_independence.cc
    test_ewald_forces_fd.cc
    test_ewald_madelung.cc
    test_ewald_spce_water.cc
    test_gcmc_detailed_balance.cc
    test_gcmc_ch4_irmof1_isotherm.cc
    test_gcmc_widom_henry.cc
    test_lattice_minimum_image.cc
    test_peng_robinson.cc
    test_philox_random123_vectors.cc
    test_lennard_jones_forces_fd.cc
    test_lennard_jones_nist.cc
    test_lennard_jones_tail_correction.cc
    test_neighbor_list_brute_force.cc
)
aleator_configure_target(aleator_validation_tests)

# Real NIST SRSW reference configurations (see data/nist_lj/PROVENANCE.md),
# read from disk rather than embedded as string literals — the largest is
# 800 particles / ~65 KB, unwieldy as a raw string. ALEATOR_TEST_DATA_DIR is
# the tests/validation/ source directory; test_lennard_jones_nist.cc joins
# "data/nist_lj/..." onto it.
target_compile_definitions(aleator_validation_tests
    PRIVATE
        ALEATOR_TEST_DATA_DIR="${CMAKE_CURRENT_SOURCE_DIR}"
)

target_link_libraries(aleator_validation_tests
    PRIVATE
        Catch2::Catch2WithMain
        aleator::core::common
        aleator::core::geometry
        aleator::core::neighbor
        aleator::core::memory
        aleator::core::math
        aleator::engines::monte_carlo
        aleator::forcefield::electrostatics
        aleator::forcefield::pairwise
        aleator::io
)

include(Catch)
catch_discover_tests(aleator_validation_tests PROPERTIES LABELS "validation" DISCOVERY_MODE PRE_TEST)
