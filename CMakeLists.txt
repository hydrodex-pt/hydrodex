cmake_minimum_required(VERSION 3.24...3.31)

project(
    hydrodex 
    VERSION
        1.0.0
    DESCRIPTION
        "Compute fluid profiles and the gravitational wave spectrum for direct and inverse phase transitions."
    LANGUAGES 
        CXX
)

# use Release as default build type
if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
    message(STATUS "Setting build type to 'Release' as none was specified.")
    set(CMAKE_BUILD_TYPE "Release" CACHE STRING "Choose the type of build." FORCE)
    set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS "Debug" "Release" "MinSizeRel" "RelWithDebInfo")
endif()

# Build config
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_MACOSX_RPATH ON)

# hydrodex build options
option(BUILD_PYTHONLIB "Build the Python interface" ON)
option(BUILD_DOCS "Build Doxygen documentation" ON)
option(BUILD_EXECUTABLE "Build the hydrodex executable" ON)

# cmake compiler options
option(CMAKE_POSITION_INDEPENDENT_CODE "Enable PIC" ON) # -fPIC
option(CMAKE_INTERPROCEDURAL_OPTIMIZATION "Enable IPO/LTO" ON) # -flto
# compiler options
option(OPTIMIZE_NATIVE "Optimize code for the machine architecture (-march=native)" ON)
option(DEVELOPMENT_WARNINGS "Switch on warnings for development (-Wall)" OFF)
option(VECTORIZATION_REPORTS "Enable loop vectorization reports (-Rpass)" OFF)
# option for reducing doxygen output
option(QUIET_DOXYGEN "Only display warnings and errors when creating the documentation" ON)

# Required dependencies
find_package(Boost CONFIG REQUIRED)
find_package(GSL REQUIRED)
find_package(OpenMP REQUIRED)
# Optional dependencies
find_package(Python3 COMPONENTS Interpreter Development.Module)
find_package(pybind11)
if(NOT (Python3_FOUND AND pybind11_FOUND))
    if(BUILD_PYTHONLIB)
        message(FATAL_ERROR
            "Python3 and/or pybind11 not found."
            " Please either provide help for locating the packages (e.g."
            " Python3_ROOT_DIR and pybind11_DIR) or set BUILD_PYTHONLIB=OFF."
        )
    else()
        message(STATUS
            "Python3 and/or pybind11 not found."
            " You will not be able to build the Python interface."
        )
    endif()
endif()
find_package(Doxygen 1.9.7)
if(NOT Doxygen_FOUND)
    if(BUILD_DOCS)
        message(FATAL_ERROR
            "Doxygen (>= 1.9.7) not found."
            " Please either provide help for locating the package (e.g."
            " DOXYGEN_EXECUTABLE) or set BUILD_DOCS=OFF."
        )
    else()
        message(STATUS
            "Doxygen (>= 1.9.7) not found. You will not be able to build the"
            " documentation."
        )
    endif()
endif()

# set compiler options
add_library(hydrodex_compiler_flags INTERFACE)
if(CMAKE_CXX_COMPILER_ID MATCHES "Clang|GNU")
    if (OPTIMIZE_NATIVE)
        message(STATUS "Code is optimized for the machine-specific architecture.")
        target_compile_options(hydrodex_compiler_flags INTERFACE -march=native)
    endif()
    if (DEVELOPMENT_WARNINGS)
        message(STATUS "Compilation warnings are enabled.")
        target_compile_options(hydrodex_compiler_flags INTERFACE -Wall)
    else()
        target_compile_options(hydrodex_compiler_flags INTERFACE -Wno-abi)
        target_link_options(hydrodex_compiler_flags INTERFACE -Wno-abi)
    endif()
    if (VECTORIZATION_REPORTS)
        message(STATUS "Vectorization reports are enabled.")
        target_compile_options(hydrodex_compiler_flags INTERFACE
            -Rpass=loop-vectorize
            -Rpass-missed=loop-vectorize
        )
    endif()
else()
    message(
        WARNING 
        "Compiler '${CMAKE_CXX_COMPILER_ID}' is not explicitly supported (supported compilers are GNU and Clang).\n"
        "The compilation flag options DEVELOPMENT_WARNINGS, OPTIMIZE_NATIVE, and VECTORIZATION_REPORTS are ignored.\n"
        "If you want to enable similar flags, please pass them manually using:\n"
        "    -DCMAKE_CXX_FLAGS=\"<your flags>\"\n"
        "or set environment variables like:\n"
        "    CXXFLAGS=\"<your flags>\"\n"
    )
endif()

# Output directories
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)
# Python site-packages directory (relative to lib)
set(PYSITE_SUBDIRECTORY python${Python3_VERSION_MAJOR}.${Python3_VERSION_MINOR}/site-packages)

# define variable for setting relative RPATH
if(APPLE)
    set(RPATH_EXECUTABLE "@executable_path")
    set(RPATH_LIBRARY "@loader_path")
elseif(UNIX)
    set(RPATH_EXECUTABLE "\\$ORIGIN")
    set(RPATH_LIBRARY "\\$ORIGIN")
endif()

# Generate version header from the project version
if(PROJECT_VERSION_PATCH STREQUAL "")
    set(HYDRODEX_VERSION_PATCH 0)
else()
    set(HYDRODEX_VERSION_PATCH ${PROJECT_VERSION_PATCH})
endif()
set(GENERATED_INCLUDE_DIR ${CMAKE_BINARY_DIR}/include)
configure_file(
    ${PROJECT_SOURCE_DIR}/include/hydrodex/version.hpp.in
    ${GENERATED_INCLUDE_DIR}/hydrodex/version.hpp
    @ONLY
)

# C++ library (linked statically into the Python module when installing via pip)
if(SKBUILD)
    set(HYDRODEX_LIBRARY_TYPE STATIC)
else()
    set(HYDRODEX_LIBRARY_TYPE SHARED)
endif()
add_library(hydrodex ${HYDRODEX_LIBRARY_TYPE}
    src/fluid_profile.cpp
    src/profile_calculator.cpp
    src/sound_shell_model.cpp
    src/sound_shell_spectrum.cpp
    src/utils.cpp
    src/settings.cpp
    src/version.cpp
)
target_include_directories(hydrodex
    PUBLIC ${PROJECT_SOURCE_DIR}/include
    PUBLIC ${GENERATED_INCLUDE_DIR}
)
target_link_libraries(hydrodex
    PUBLIC Boost::headers
    PUBLIC GSL::gsl GSL::gslcblas
    PUBLIC OpenMP::OpenMP_CXX
    PRIVATE hydrodex_compiler_flags
)
# Install shared library + headers (unless building via pip)
if(NOT SKBUILD)
    install(TARGETS hydrodex
        LIBRARY DESTINATION lib
        ARCHIVE DESTINATION lib
    )
    install(DIRECTORY ${PROJECT_SOURCE_DIR}/include/ DESTINATION include
        PATTERN "*.in" EXCLUDE
    )
    install(DIRECTORY ${GENERATED_INCLUDE_DIR}/ DESTINATION include)
    # license, readme and changelog
    install(FILES
        ${PROJECT_SOURCE_DIR}/LICENSE
        ${PROJECT_SOURCE_DIR}/README.md
        ${PROJECT_SOURCE_DIR}/CHANGELOG.md
        DESTINATION share/${PROJECT_NAME}/docs
    )
endif()

# Python interface
if(Python3_FOUND AND pybind11_FOUND)
    pybind11_add_module(pythonlib
        src/pythonlib/fluid_profile_python.cpp
        src/pythonlib/profile_calculator_python.cpp
        src/pythonlib/sound_shell_model_python.cpp
        src/pythonlib/sound_shell_spectrum_python.cpp
        src/pythonlib/utils_python.cpp
        src/pythonlib/settings_python.cpp
        src/pythonlib/hydrodex.cpp
    )
    target_include_directories(pythonlib 
        PRIVATE ${PROJECT_SOURCE_DIR}/src/pythonlib
    )
    target_link_libraries(pythonlib 
        PRIVATE hydrodex hydrodex_compiler_flags
    )
    set_target_properties(pythonlib PROPERTIES
        OUTPUT_NAME "hydrodex"
        LIBRARY_OUTPUT_DIRECTORY ${CMAKE_LIBRARY_OUTPUT_DIRECTORY}/${PYSITE_SUBDIRECTORY}
    )
    if(NOT SKBUILD)
        set_target_properties(pythonlib PROPERTIES
            BUILD_RPATH "${CMAKE_LIBRARY_OUTPUT_DIRECTORY};${RPATH_EXECUTABLE}/../.."
            INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib;${RPATH_EXECUTABLE}/../.."
        )
    endif()
    message(STATUS "BUILD_PYTHONLIB for building the Python interface is set to ${BUILD_PYTHONLIB}.")
    if(SKBUILD)
        # install into the wheel when installing via pip
        # (no RPATH needed since libhydrodex is linked statically)
        message(STATUS "The Python interface will be installed into ${SKBUILD_PLATLIB_DIR}.")
        set_target_properties(pythonlib PROPERTIES
            BUILD_WITH_INSTALL_RPATH TRUE
        )
        install(TARGETS pythonlib
            LIBRARY DESTINATION .
            COMPONENT python
        )
    elseif(BUILD_PYTHONLIB)
        # install
        message(STATUS "The Python interface will be built and installed to ${CMAKE_INSTALL_PREFIX}/lib/${PYSITE_SUBDIRECTORY}.")
        install(TARGETS pythonlib
            LIBRARY DESTINATION lib/${PYSITE_SUBDIRECTORY}
        )
    else()
        # do not build pythonlib by default unless requested
        message(STATUS "The Python interface will not be built by default and will not be installed.")
        set_property(TARGET pythonlib PROPERTY EXCLUDE_FROM_ALL TRUE)
    endif()
else()
    message(STATUS "Python or pybind11 were not found. The Python interface cannot be built.")
endif()

# Executable (target name differs from the library target, output name is hydrodex)
add_executable(hydrodex_bin cli/hydrodex.cpp)
target_link_libraries(hydrodex_bin
    PRIVATE hydrodex hydrodex_compiler_flags
)
set_target_properties(hydrodex_bin PROPERTIES
    OUTPUT_NAME "hydrodex"
    BUILD_RPATH "${CMAKE_LIBRARY_OUTPUT_DIRECTORY};${RPATH_EXECUTABLE}/../lib"
    INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib;${RPATH_EXECUTABLE}/../../lib"
)
message(STATUS "BUILD_EXECUTABLE for building the hydrodex executable is set to ${BUILD_EXECUTABLE}")
if(BUILD_EXECUTABLE)
    message(STATUS "The hydrodex executable will be built and installed to ${CMAKE_INSTALL_PREFIX}/bin/hydrodex.")
    install(TARGETS hydrodex_bin DESTINATION bin/hydrodex)
else()
    # do not build executable by default unless requested
    message(STATUS "The hydrodex executable will not be built by default and will not be installed.")
    set_property(TARGET hydrodex_bin PROPERTY EXCLUDE_FROM_ALL TRUE)
endif()

# Documentation
if(DOXYGEN_FOUND)
    # prepare config file
    set(DOXYGEN_CONFIG_OUT ${CMAKE_BINARY_DIR}/doxygen.config)
    file(READ ${PROJECT_SOURCE_DIR}/docs/doxygen.config DOXYGEN_CONFIG_CONTENTS)
    file(WRITE ${DOXYGEN_CONFIG_OUT} "${DOXYGEN_CONFIG_CONTENTS}\n")
    file(APPEND ${DOXYGEN_CONFIG_OUT} "OUTPUT_DIRECTORY = ${CMAKE_BINARY_DIR}/docs\n")
    file(APPEND ${DOXYGEN_CONFIG_OUT} "INPUT = ${PROJECT_SOURCE_DIR}/README.md ${PROJECT_SOURCE_DIR}/CHANGELOG.md ${PROJECT_SOURCE_DIR}/include/hydrodex ${GENERATED_INCLUDE_DIR}/hydrodex ${PROJECT_SOURCE_DIR}/src\n")
    file(APPEND ${DOXYGEN_CONFIG_OUT} "USE_MDFILE_AS_MAINPAGE = ${PROJECT_SOURCE_DIR}/README.md\n")
    if(QUIET_DOXYGEN)
        file(APPEND ${DOXYGEN_CONFIG_OUT} "QUIET = YES\n")
    endif()
    # Doxygen docs target
    add_custom_target(docs
        COMMAND  Doxygen::doxygen ${DOXYGEN_CONFIG_OUT}
        # images referenced in the README via <img> are not copied by Doxygen
        COMMAND  ${CMAKE_COMMAND} -E copy_directory ${PROJECT_SOURCE_DIR}/docs/images ${CMAKE_BINARY_DIR}/docs/html/docs/images
        WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
        COMMENT "Generating API documentation into build/docs"
    )
    # Sphinx documentation of the Python interface (read from the compiled module)
    if(BUILD_PYTHONLIB AND TARGET pythonlib AND NOT SKBUILD)
        execute_process(
            COMMAND ${Python3_EXECUTABLE} -c "import sphinx, furo"
            RESULT_VARIABLE SPHINX_IMPORT_RESULT
            OUTPUT_QUIET ERROR_QUIET
        )
        if(SPHINX_IMPORT_RESULT EQUAL 0)
            message(STATUS "The documentation will include the Python interface.")
            add_custom_target(docs_python
                COMMAND ${CMAKE_COMMAND} -E env HYDRODEX_MODULE=$<TARGET_FILE:pythonlib>
                    ${Python3_EXECUTABLE} -m sphinx -q -E -b html
                    -d ${CMAKE_BINARY_DIR}/sphinx_doctrees
                    ${PROJECT_SOURCE_DIR}/docs/sphinx ${CMAKE_BINARY_DIR}/docs/html/python
                WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
                COMMENT "Generating Python API documentation into build/docs/html/python"
            )
            add_dependencies(docs_python pythonlib)
            add_dependencies(docs docs_python)
            # layout with an additional navigation tab linking to the Python docs;
            # start from the full default layout of the installed Doxygen version,
            # since Doxygen leaves out all parts missing from a custom layout file
            set(DOXYGEN_LAYOUT_OUT ${CMAKE_BINARY_DIR}/DoxygenLayout.xml)
            execute_process(
                COMMAND ${DOXYGEN_EXECUTABLE} -l ${DOXYGEN_LAYOUT_OUT}
                WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
                OUTPUT_QUIET
            )
            file(READ ${DOXYGEN_LAYOUT_OUT} DOXYGEN_LAYOUT_CONTENTS)
            string(REPLACE "</navindex>"
                "  <tab type=\"user\" url=\"python/index.html\" title=\"Python API\"/>\n  </navindex>"
                DOXYGEN_LAYOUT_CONTENTS "${DOXYGEN_LAYOUT_CONTENTS}"
            )
            file(WRITE ${DOXYGEN_LAYOUT_OUT} "${DOXYGEN_LAYOUT_CONTENTS}")
            file(APPEND ${DOXYGEN_CONFIG_OUT} "LAYOUT_FILE = ${DOXYGEN_LAYOUT_OUT}\n")
        else()
            message(STATUS
                "Sphinx and/or furo not found for ${Python3_EXECUTABLE}."
                " The documentation will not include the Python interface."
            )
        endif()
    endif()
    message(STATUS "BUILD_DOCS for building the documentation is set to ${BUILD_DOCS}")
    # build docs by default (and install them) when requested
    if(BUILD_DOCS)
        message(STATUS "The documentation will be built and installed to ${CMAKE_INSTALL_PREFIX}/share/${PROJECT_NAME}/.")
        set_property(TARGET docs PROPERTY EXCLUDE_FROM_ALL FALSE)
        install(DIRECTORY ${CMAKE_BINARY_DIR}/docs DESTINATION share/${PROJECT_NAME}/)
    else()
        message(STATUS "The documentation will not be built by default and will not be installed.")
    endif()
endif()

# Example Jupyter (Python) notebooks
if(BUILD_PYTHONLIB AND NOT SKBUILD)
    install(DIRECTORY ${PROJECT_SOURCE_DIR}/examples/
        DESTINATION share/${PROJECT_NAME}/examples
        FILES_MATCHING PATTERN "*.ipynb"
    )
endif()

# target to build everything
add_custom_target(everything DEPENDS hydrodex hydrodex_bin)
if(Python3_FOUND AND pybind11_FOUND)
    add_dependencies(everything pythonlib)
endif()
if(DOXYGEN_FOUND)
    add_dependencies(everything docs)
endif()
