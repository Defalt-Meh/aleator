add_library(aleator_core_geometry STATIC
    lattice.cc
)
add_library(aleator::core::geometry ALIAS aleator_core_geometry)
target_sources(aleator_core_geometry PUBLIC
    FILE_SET HEADERS
    BASE_DIRS ${PROJECT_SOURCE_DIR}/src
    FILES
        lattice.hpp
)
aleator_configure_target(aleator_core_geometry)

target_link_libraries(aleator_core_geometry PRIVATE aleator::core::common)
