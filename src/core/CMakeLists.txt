# core/ depends on nothing above it (CLAUDE.md #3) — no target in this tree
# may link against forcefield::*, engines::*, io::*, or cli::*.

add_library(aleator_core_common INTERFACE)
add_library(aleator::core::common ALIAS aleator_core_common)
target_sources(aleator_core_common INTERFACE
    FILE_SET HEADERS
    BASE_DIRS ${PROJECT_SOURCE_DIR}/src
    FILES exceptions.hpp
)
aleator_configure_interface_target(aleator_core_common)

add_subdirectory(memory)
add_subdirectory(math)
add_subdirectory(geometry)
add_subdirectory(neighbor)
