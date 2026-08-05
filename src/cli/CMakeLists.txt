# cli/ is the top of the dependency graph (CLAUDE.md #3): the only target
# permitted to link core/, forcefield/, engines/, and io/ together.

add_executable(aleator_cli main.cc)
set_target_properties(aleator_cli PROPERTIES OUTPUT_NAME aleator)
aleator_configure_target(aleator_cli)

target_link_libraries(aleator_cli
    PRIVATE
        aleator::version
        aleator::core::common
        aleator::core::geometry
        aleator::core::math
        aleator::core::neighbor
        aleator::forcefield::pairwise
        aleator::engines::monte_carlo
        aleator::io
        spdlog::spdlog
)

install(TARGETS aleator_cli RUNTIME DESTINATION bin)
