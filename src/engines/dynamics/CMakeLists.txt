add_library(aleator_engines_dynamics STATIC
    integrator.cc
)
add_library(aleator::engines::dynamics ALIAS aleator_engines_dynamics)
target_sources(aleator_engines_dynamics PUBLIC
    FILE_SET HEADERS
    BASE_DIRS ${PROJECT_SOURCE_DIR}/src
    FILES integrator.hpp
)
aleator_configure_target(aleator_engines_dynamics)

target_link_libraries(aleator_engines_dynamics
    PUBLIC
        aleator::core::geometry
        aleator::core::math
        aleator::forcefield::common
    PRIVATE
        aleator::core::common
)
