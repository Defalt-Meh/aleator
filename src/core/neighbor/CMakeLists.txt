add_library(aleator_core_neighbor STATIC
    cell_list.cc
    verlet_list.cc
)
add_library(aleator::core::neighbor ALIAS aleator_core_neighbor)
target_sources(aleator_core_neighbor PUBLIC
    FILE_SET HEADERS
    BASE_DIRS ${PROJECT_SOURCE_DIR}/src
    FILES
        cell_list.hpp
        neighbor_list.hpp
        verlet_list.hpp
)
aleator_configure_target(aleator_core_neighbor)

target_link_libraries(aleator_core_neighbor
    PUBLIC
        aleator::core::geometry
        aleator::core::math
)
