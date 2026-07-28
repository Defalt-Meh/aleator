add_library(aleator_engines_monte_carlo STATIC
    monte_carlo_engine.cc
    peng_robinson.cc
)
add_library(aleator::engines::monte_carlo ALIAS aleator_engines_monte_carlo)
target_sources(aleator_engines_monte_carlo PUBLIC
    FILE_SET HEADERS
    BASE_DIRS ${PROJECT_SOURCE_DIR}/src
    FILES monte_carlo_engine.hpp peng_robinson.hpp molecule_species.hpp gcmc_acceptance.hpp
)
aleator_configure_target(aleator_engines_monte_carlo)

target_link_libraries(aleator_engines_monte_carlo
    PUBLIC
        aleator::core::geometry
        aleator::core::math
        aleator::forcefield::common
    PRIVATE
        aleator::core::common
)
