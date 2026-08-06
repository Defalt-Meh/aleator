add_library(aleator_engines_geometry_analysis STATIC
    pore_analysis.cc
    periodic_graph.cc
    voronoi_network.cc
)
add_library(aleator::engines::geometry_analysis ALIAS aleator_engines_geometry_analysis)
target_sources(aleator_engines_geometry_analysis PUBLIC
    FILE_SET HEADERS
    BASE_DIRS ${PROJECT_SOURCE_DIR}/src
    FILES pore_analysis.hpp periodic_graph.hpp voronoi_network.hpp
)
aleator_configure_target(aleator_engines_geometry_analysis)

target_link_libraries(aleator_engines_geometry_analysis
    PUBLIC
        aleator::core::geometry
        aleator::core::math
    PRIVATE
        aleator::core::common
        VORO::voro++
)
