add_library(aleator_core_geometry STATIC
    lattice.cc
    supercell.cc
)
add_library(aleator::core::geometry ALIAS aleator_core_geometry)
target_sources(aleator_core_geometry PUBLIC
    FILE_SET HEADERS
    BASE_DIRS ${PROJECT_SOURCE_DIR}/src
    FILES
        lattice.hpp
        supercell.hpp
)
aleator_configure_target(aleator_core_geometry)

# supercell.hpp is a public header that includes core/math/particle_data.hpp,
# so that dependency must be PUBLIC (same rationale as core/math's own
# PUBLIC dependency on core/memory, just one layer up).
target_link_libraries(aleator_core_geometry
    PUBLIC
        aleator::core::math
    PRIVATE
        aleator::core::common
)
