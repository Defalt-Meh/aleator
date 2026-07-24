# Warnings-as-errors interface target, linked PRIVATE into every Aleator
# target. Kept in one place so the warning set is identical across GCC,
# Clang, and MSVC rather than drifting per-target.
add_library(aleator_warnings INTERFACE)
add_library(aleator::warnings ALIAS aleator_warnings)

if(CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
    target_compile_options(aleator_warnings INTERFACE
        -Wall
        -Wextra
        -Wpedantic
        -Wshadow
        -Wconversion
        -Wsign-conversion
        -Wcast-qual
        -Wold-style-cast
        -Wnon-virtual-dtor
        -Woverloaded-virtual
        -Wdouble-promotion
        -Werror
    )
elseif(CMAKE_CXX_COMPILER_ID STREQUAL "MSVC")
    target_compile_options(aleator_warnings INTERFACE
        /W4
        /permissive-
        /WX
    )
endif()
