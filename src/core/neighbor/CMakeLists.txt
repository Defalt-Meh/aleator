# Pure abstract interface: no translation unit.
add_library(aleator_core_neighbor INTERFACE)
add_library(aleator::core::neighbor ALIAS aleator_core_neighbor)
target_sources(aleator_core_neighbor INTERFACE
    FILE_SET HEADERS
    BASE_DIRS ${PROJECT_SOURCE_DIR}/src
    FILES
        neighbor_list.hpp
)
aleator_configure_interface_target(aleator_core_neighbor)

target_link_libraries(aleator_core_neighbor
    INTERFACE
        aleator::core::geometry
        aleator::core::math
)
