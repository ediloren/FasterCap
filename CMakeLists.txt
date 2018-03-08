# minimum version is 3.0 to be able to handle Eigen library
cmake_minimum_required (VERSION 2.8)
project (FasterCap)

# main project specific settings
#
# option to compile as headless, or with GUI
set(FASTFIELDSOLVERS_HEADLESS OFF CACHE BOOL "Compile as headless (i.e. no GUI)")

#
# Global Compiler and Linker Settings
#

#
# using only GCC, therefore configuring build only for GCC
#
# .cmake files are script files; not using here any script
#include(cMake/custom_include.cmake)
#configure_file(config.h.cmake ${CMAKE_CURRENT_BINARY_DIR}/config.h)
#
# if headless, pass on the "#define FCG_HEADLESS"
if(${FASTFIELDSOLVERS_HEADLESS} MATCHES ON)
    add_definitions("-DFCG_HEADLESS -DwxUSE_GUI=0")
else()
    add_definitions("-DwxUSE_GUI=1")
endif()
#
# pass to gcc the argument -DHAVE_CONFIG_H . This defines the pre-processor token HAVE_CONFIG_H
# exactly as if you had #define HAVE_CONFIG_H right at the start of each of your source files.
# HAVE_CONFIG_H is a pre-defined macro in the build-system that generates a config.h
# so that code knows whether it needs to #include config.h
#add_definitions(-DHAVE_CONFIG_H)
#
# List of options to pass to the compiler (note this is a CMake list. Do not enclose in ")
set(FASTFIELDSOLVERS_COMPILE_OPTIONS ${FASTFIELDSOLVERS_COMPILE_OPTIONS} -Wextra -Wall -fopenmp)
# -fopenmp: use OpenMP pragmas
# -lm: link the math library
# -lstdc++: link the standard c++ library (needed for headless)
set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -fopenmp -lm -lstdc++")
# -s: strip all symbols from the executable
set(CMAKE_EXE_LINKER_FLAGS_RELEASE "${CMAKE_EXE_LINKER_FLAGS_RELEASE} -s")

# wxWidgets usage
# See CMake docs "FindwxWidgets"
# In general, for *nix CMake already generates the 'wx-confg' options that are needed
# for compiling and linking with wxWidgets. However, if anything special is needed,
# it must be tweaked into "wxWidgets_CONFIG_OPTIONS".
# In this case, we pass the options `wx-config --version=3.0 --static=no --debug`,
# but '--debug' is only needed in the debug build. Therefore, we use a CMake generator-expression
# to define the flag "--debug" only in case of debug build (this is nested: the first clause
# $<CONFIG:Debug> generates 1 or 0 if CONFIG is Debug or not. Then we use the conditional
# expression $<condition:true_string>, so if the previous evaluates to 1, here the string is "--debug",
# otherwise nothing.
set(wxWidgets_CONFIG_OPTIONS --version=3.0 --static=no $<$<CONFIG:Debug>:--debug>)
# these flags are BUILD_TYPE dependent
# -D__WXDEBUG__: define that wxWidgets are in debug version
#set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -D__WXDEBUG__")
if(${FASTFIELDSOLVERS_HEADLESS} MATCHES ON)
    # we only need the 'base' wxWidget package for headless compile
    find_package(wxWidgets COMPONENTS base)
else()
    # probably not all required. Note however that order is important.
    #find_package(wxWidgets COMPONENTS html xrc net core base msw png jpeg tiff zlib kernel32 user32 gdi32 winspool comdlg32 advapi32 shell32 ole32 oleaut32 uuid comctl32 sock32 odbc32 iphlpapi)
    find_package(wxWidgets COMPONENTS core base html)
endif()
# Configuration to use
# (e.g., msw, mswd, mswu, mswunivud, etc.)
# wxWidgets_CONFIGURATION
# general include (needed for wxWidgets in CMake), see CMake docs "FindwxWidgets"
include(${wxWidgets_USE_FILE})
# add to WXWIDGETS_LIBS the wxWidgets libraries
# we'll subsequently link WXWIDGETS_LIBS with 'target_link_libraries'
set(WXWIDGETS_LIBS ${WXWIDGETS_LIBS} ${wxWidgets_LIBRARIES})

# general include_directories, as LinAlgebra and Geometry
# are one level under the base path, and they are included
# in the sources as "#include "LinAlgebra/Mtx.h" etc.
include_directories("..")
# same as above for the main includes in the base path; could use ${CMAKE_SOURCE_DIR}
include_directories(".")

include_directories ("../LinAlgebra")
add_subdirectory("../LinAlgebra" LinAlgebra)
include_directories ("../Geometry")
add_subdirectory("../Geometry" Geometry)
set (EXTRA_LIBS ${EXTRA_LIBS} LinAlgebra Geometry)

include_directories ("${PROJECT_SOURCE_DIR}/Solver")
include("./Solver/CMakeLists.txt")
# do not create a library, to speed up compilation in case of file modifications
# (only obj to recompile and final linking, no need to re-create also the library)
#add_subdirectory(Solver)
#set (EXTRA_LIBS ${EXTRA_LIBS} Solver)

set(FastFieldSolvers_Sources ${FastFieldSolvers_Sources}
    FasterCapConsole.cpp
    FasterCapConsole.h
    FasterCapGlobal.cpp
    FasterCapGlobal.h
    test.cpp
    test.h)

# if running with a GUI
if(${FASTFIELDSOLVERS_HEADLESS} MATCHES OFF)
    set(FastFieldSolvers_Sources ${FastFieldSolvers_Sources}
        AboutBox.cpp
        AboutBox.h
        FasterCapApp.cpp
        FasterCapApp.h
        FasterCapMain.cpp
        FasterCapMain.h
        LicenseBox.cpp
        LicenseBox.h
        RunDialog.cpp
        RunDialog.h
        resource.rc)
    # and if also running on Windows
    if(CMAKE_SYSTEM_NAME STREQUAL Windows)
        set(FastFieldSolvers_Sources ${FastFieldSolvers_Sources}
            AutomationHelper.cpp
            AutomationHelper.h)
    endif()
endif()

add_executable (FasterCap ${FastFieldSolvers_Sources} ${FasterCap_Sources})
# Specify compile options to use when compiling a given target.
target_compile_options(FasterCap PUBLIC ${FASTFIELDSOLVERS_COMPILE_OPTIONS})
target_link_libraries(FasterCap ${EXTRA_LIBS} ${WXWIDGETS_LIBS} ${EIGEN_LIBS})
