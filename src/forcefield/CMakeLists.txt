# forcefield/ depends only on core/ (CLAUDE.md #3) — nothing here may link
# against engines::*, io::*, or cli::*.

add_library(aleator_forcefield_common INTERFACE)
add_library(aleator::forcefield::common ALIAS aleator_forcefield_common)
target_sources(aleator_forcefield_common INTERFACE
    FILE_SET HEADERS
    BASE_DIRS ${PROJECT_SOURCE_DIR}/src
    FILES force_field.hpp
)
aleator_configure_interface_target(aleator_forcefield_common)
target_link_libraries(aleator_forcefield_common
    INTERFACE
        aleator::core::geometry
        aleator::core::math
        aleator::core::memory
        aleator::core::neighbor
)

add_subdirectory(pairwise)
add_subdirectory(electrostatics)
add_subdirectory(parameters)
