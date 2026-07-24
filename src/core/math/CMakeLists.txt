add_library(aleator_core_math STATIC
    simd_kernels.cc
    counter_based_rng.cc
)
add_library(aleator::core::math ALIAS aleator_core_math)
target_sources(aleator_core_math PUBLIC
    FILE_SET HEADERS
    BASE_DIRS ${PROJECT_SOURCE_DIR}/src
    FILES
        simd_kernels.hpp
        particle_data.hpp
        counter_based_rng.hpp
)
aleator_configure_target(aleator_core_math)

# particle_data.hpp is a public header that includes
# core/memory/aligned_allocator.hpp, so that dependency must be PUBLIC.
target_link_libraries(aleator_core_math
    PUBLIC
        aleator::core::memory
    PRIVATE
        aleator::core::common
        hwy::hwy
)
