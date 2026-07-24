add_library(aleator_forcefield_parameters STATIC
    parameter_loader.cc
)
add_library(aleator::forcefield::parameters ALIAS aleator_forcefield_parameters)
target_sources(aleator_forcefield_parameters PUBLIC
    FILE_SET HEADERS
    BASE_DIRS ${PROJECT_SOURCE_DIR}/src
    FILES parameter_loader.hpp
)
aleator_configure_target(aleator_forcefield_parameters)

target_link_libraries(aleator_forcefield_parameters
    PUBLIC
        aleator::forcefield::pairwise
    PRIVATE
        aleator::core::common
)
