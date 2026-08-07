# engines/ depend on core/ + forcefield/ (CLAUDE.md #3). No lateral includes
# between engines: monte_carlo and geometry_analysis must not link against
# each other.

add_subdirectory(monte_carlo)
add_subdirectory(geometry_analysis)
