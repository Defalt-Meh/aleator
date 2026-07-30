add_executable(aleator_integration_tests
    test_cli_end_to_end.cc
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

# test_cli_end_to_end.cc drives the real `aleator` binary as a subprocess --
# it needs that binary actually built first, and needs to know where it
# landed and where examples/ lives (both vary by platform/build config).
add_dependencies(aleator_integration_tests aleator_cli)
target_compile_definitions(aleator_integration_tests
    PRIVATE
        ALEATOR_CLI_BINARY_PATH="$<TARGET_FILE:aleator_cli>"
        ALEATOR_EXAMPLES_DIR="${PROJECT_SOURCE_DIR}/examples"
)

include(Catch)
catch_discover_tests(aleator_integration_tests PROPERTIES LABELS "integration" DISCOVERY_MODE PRE_TEST)
