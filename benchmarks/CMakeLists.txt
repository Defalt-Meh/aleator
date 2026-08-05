add_executable(aleator_bench
    bench_neighbor_list.cc
    bench_vector_sum.cc
    bench_gcmc_ch4_irmof1.cc
)
aleator_configure_target(aleator_bench)

target_link_libraries(aleator_bench
    PRIVATE
        benchmark::benchmark
        benchmark::benchmark_main
        aleator::core::math
        aleator::core::neighbor
        aleator::engines::monte_carlo
        aleator::forcefield::pairwise
)

# ALEATOR_BENCH_NATIVE scopes -march=native to this executable only (set by
# the "bench" preset). It must never reach a target that ships in a wheel or
# release binary — see CLAUDE.md #2.3 and the top-level CMakeLists.txt guard.
if(ALEATOR_BENCH_NATIVE)
    if(CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
        target_compile_options(aleator_bench PRIVATE -march=native)
    else()
        message(WARNING "ALEATOR_BENCH_NATIVE has no effect for ${CMAKE_CXX_COMPILER_ID} (GCC/Clang only).")
    endif()
endif()
