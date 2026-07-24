# bindings/ is the Python-facing edge of the dependency graph: it may link
# core/ (and, once real physics exists, forcefield/engines/io) but nothing
# links back to it. Only built when ALEATOR_BUILD_PYTHON is ON.

if(CMAKE_VERSION VERSION_LESS 3.18)
    set(ALEATOR_PYTHON_DEV_COMPONENT Development)
else()
    set(ALEATOR_PYTHON_DEV_COMPONENT Development.Module)
endif()

find_package(Python 3.9 COMPONENTS Interpreter ${ALEATOR_PYTHON_DEV_COMPONENT} REQUIRED)

if(DEFINED SKBUILD)
    # scikit-build-core wheel/editable build: fetch nanobind from PyPI (see
    # [build-system] requires in pyproject.toml) rather than from vcpkg.
    # vcpkg's "python" manifest feature is disabled for this build
    # (-DVCPKG_MANIFEST_NO_DEFAULT_FEATURES=ON, also set in pyproject.toml)
    # precisely so its bundled Python 3 never gets installed and never gets
    # a chance to hijack the find_package(Python) call above — the
    # extension must match the exact host interpreter's ABI, which is what
    # PyPI's nanobind (headers usable with any host Python) guarantees.
    execute_process(
        COMMAND "${Python_EXECUTABLE}" -m nanobind --cmake_dir
        OUTPUT_STRIP_TRAILING_WHITESPACE
        OUTPUT_VARIABLE nanobind_ROOT
    )
endif()
find_package(nanobind CONFIG REQUIRED)

nanobind_add_module(_aleator
    STABLE_ABI
    module.cc
)
target_compile_features(_aleator PRIVATE cxx_std_23)
set_target_properties(_aleator PROPERTIES CXX_EXTENSIONS OFF)

# Deliberately NOT linking aleator::warnings here: module.cc pulls in
# CPython's own headers (via nanobind), which are riddled with old-style
# C casts we don't control. aleator::warnings' -Werror would be fighting
# a third-party API surface, not our code.
target_link_libraries(_aleator
    PRIVATE
        aleator::version
        aleator::core::math
)

# scikit-build-core sets CMAKE_INSTALL_PREFIX to the wheel staging root; the
# extension must land next to python/aleator/__init__.py. For a standalone
# (non-wheel) CMake build, this installs relative to CMAKE_INSTALL_PREFIX and
# is otherwise inert.
install(TARGETS _aleator LIBRARY DESTINATION aleator)
