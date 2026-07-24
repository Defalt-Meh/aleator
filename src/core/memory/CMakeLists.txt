# Header-only: aligned allocation is a template + a small inline class, no
# separate translation unit.
add_library(aleator_core_memory INTERFACE)
add_library(aleator::core::memory ALIAS aleator_core_memory)
target_sources(aleator_core_memory INTERFACE
    FILE_SET HEADERS
    BASE_DIRS ${PROJECT_SOURCE_DIR}/src
    FILES
        aligned_allocator.hpp
        arena.hpp
)
aleator_configure_interface_target(aleator_core_memory)
