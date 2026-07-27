add_library(aleator_forcefield_pairwise STATIC
    lennard_jones.cc
)
add_library(aleator::forcefield::pairwise ALIAS aleator_forcefield_pairwise)
target_sources(aleator_forcefield_pairwise PUBLIC
    FILE_SET HEADERS
    BASE_DIRS ${PROJECT_SOURCE_DIR}/src
    FILES lennard_jones.hpp
)
aleator_configure_target(aleator_forcefield_pairwise)

target_link_libraries(aleator_forcefield_pairwise
    PUBLIC
        aleator::forcefield::common
)
