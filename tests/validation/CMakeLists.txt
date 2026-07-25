add_executable(aleator_validation_tests
    test_cif_real_structures.cc
    test_ewald_not_faked.cc
    test_lattice_minimum_image.cc
)
aleator_configure_target(aleator_validation_tests)

target_link_libraries(aleator_validation_tests
    PRIVATE
        Catch2::Catch2WithMain
        aleator::core::common
        aleator::core::geometry
        aleator::core::neighbor
        aleator::forcefield::electrostatics
        aleator::io
)

include(Catch)
catch_discover_tests(aleator_validation_tests PROPERTIES LABELS "validation" DISCOVERY_MODE PRE_TEST)
