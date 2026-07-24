# ALEATOR_ENABLE_ASAN / ALEATOR_ENABLE_UBSAN are set by the "asan" CMake
# preset. Kept as a separate interface target (rather than baked into
# aleator_warnings) so sanitizer flags never leak into release/bench builds.
option(ALEATOR_ENABLE_ASAN "Build with AddressSanitizer" OFF)
option(ALEATOR_ENABLE_UBSAN "Build with UndefinedBehaviorSanitizer" OFF)

add_library(aleator_sanitizers INTERFACE)
add_library(aleator::sanitizers ALIAS aleator_sanitizers)

if(ALEATOR_ENABLE_ASAN OR ALEATOR_ENABLE_UBSAN)
    if(NOT (CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang"))
        message(FATAL_ERROR "Sanitizers are only supported with GCC or Clang.")
    endif()

    set(_aleator_sanitizer_list "")
    if(ALEATOR_ENABLE_ASAN)
        list(APPEND _aleator_sanitizer_list "address")
    endif()
    if(ALEATOR_ENABLE_UBSAN)
        list(APPEND _aleator_sanitizer_list "undefined")
    endif()
    list(JOIN _aleator_sanitizer_list "," _aleator_sanitizer_flag_value)

    target_compile_options(aleator_sanitizers INTERFACE
        -fsanitize=${_aleator_sanitizer_flag_value}
        -fno-omit-frame-pointer
        -fno-sanitize-recover=all
    )
    target_link_options(aleator_sanitizers INTERFACE
        -fsanitize=${_aleator_sanitizer_flag_value}
    )
endif()
