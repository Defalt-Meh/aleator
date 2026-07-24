add_executable(aleator_unit_tests
    test_core_types.cc
    test_simd_kernels.cc
)
aleator_configure_target(aleator_unit_tests)

target_link_libraries(aleator_unit_tests
    PRIVATE
        Catch2::Catch2WithMain
        aleator::core::common
        aleator::core::memory
        aleator::core::math
        aleator::core::geometry
        hwy::hwy
)

include(Catch)
# PRE_TEST (rather than the default POST_BUILD) defers discovery to ctest
# time. POST_BUILD runs the freshly-linked binary as part of the build graph
# itself, which on Apple Silicon trips a known AddressSanitizer/SIP
# interposition issue (Catch2 #2411) under the asan preset.
catch_discover_tests(aleator_unit_tests PROPERTIES LABELS "unit" DISCOVERY_MODE PRE_TEST)
