add_library(aleator_forcefield_electrostatics STATIC
    ewald.cc
    ewald_incremental_state.cc
)
add_library(aleator::forcefield::electrostatics ALIAS aleator_forcefield_electrostatics)
target_sources(aleator_forcefield_electrostatics PUBLIC
    FILE_SET HEADERS
    BASE_DIRS ${PROJECT_SOURCE_DIR}/src
    FILES ewald.hpp ewald_incremental_state.hpp
)
aleator_configure_target(aleator_forcefield_electrostatics)

target_link_libraries(aleator_forcefield_electrostatics
    PUBLIC
        aleator::forcefield::common
)
