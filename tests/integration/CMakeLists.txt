add_executable(aleator_integration_tests
    test_engine_wiring.cc
    test_run_config_roundtrip.cc
)
aleator_configure_target(aleator_integration_tests)

target_link_libraries(aleator_integration_tests
    PRIVATE
        Catch2::Catch2WithMain
        aleator::core::common
        aleator::core::geometry
        aleator::core::math
        aleator::forcefield::pairwise
        aleator::engines::monte_carlo
        aleator::io
)

include(Catch)
catch_discover_tests(aleator_integration_tests PROPERTIES LABELS "integration" DISCOVERY_MODE PRE_TEST)
