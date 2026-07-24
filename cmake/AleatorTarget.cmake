# Common per-target setup shared by every Aleator library. Centralizing this
# means the C++ standard, include root, and warnings/sanitizer flags cannot
# drift between the ~15 leaf targets under src/.
function(aleator_configure_target target)
    target_include_directories(${target}
        PUBLIC
            $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/src>
            $<INSTALL_INTERFACE:include>
    )
    target_compile_features(${target} PUBLIC cxx_std_23)
    set_target_properties(${target} PROPERTIES
        CXX_EXTENSIONS OFF
        POSITION_INDEPENDENT_CODE ON
    )
    target_link_libraries(${target} PRIVATE aleator::warnings aleator::sanitizers)
endfunction()

# Header-only (INTERFACE) targets have no translation unit of their own to
# apply -Werror/sanitizer flags to, so they only need the include root and
# language standard propagated to whatever links them.
function(aleator_configure_interface_target target)
    target_include_directories(${target}
        INTERFACE
            $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/src>
            $<INSTALL_INTERFACE:include>
    )
    target_compile_features(${target} INTERFACE cxx_std_23)
endfunction()
