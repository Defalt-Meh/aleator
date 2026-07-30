# io/ depends downward only (CLAUDE.md #3): core/ and forcefield/'s public
# types, never engines/ or cli/.

add_library(aleator_io STATIC
    cif_document.cc
    config.cc
    element_table.cc
    structure_io.cc
    symmetry_operation.cc
)
add_library(aleator::io ALIAS aleator_io)
target_sources(aleator_io PUBLIC
    FILE_SET HEADERS
    BASE_DIRS ${PROJECT_SOURCE_DIR}/src
    FILES
        cif_document.hpp
        cif_error.hpp
        config.hpp
        config_error.hpp
        element_table.hpp
        structure_io.hpp
        symmetry_operation.hpp
)
aleator_configure_target(aleator_io)

target_link_libraries(aleator_io
    PUBLIC
        aleator::core::geometry
        aleator::core::math
    PRIVATE
        aleator::core::common
        tomlplusplus::tomlplusplus
)
