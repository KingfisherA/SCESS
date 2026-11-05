# Build script for SCESS (Schrödinger's Cat's Electronic System Solver)
#
#   cmake -S . -B build -G Ninja -DCMAKE_BUILD_TYPE=Release
#   cmake --build build
#
# Outputs, one folder per target:
#   build/QuantumSCESS   libQuantumSCESS, the shared library with the quantum chemistry code
#   build/SCESS          SCESS, the command line driver
#   build/ViewerSCESS    ViewerSCESS, the graphical viewer
#                        (macOS: ViewerSCESS.app via cmake --build build --target ViewerSCESSApp)
#   build/Dahlia         libDahlia, the graphics engine (built from SCESS_DAHLIA_DIR)

# 3.24 for $<LINK_LIBRARY:FRAMEWORK,...> (also Dahlia's minimum)
cmake_minimum_required(VERSION 3.24...4.0)

project(SCESS VERSION 1.0.0 LANGUAGES C CXX)
if(APPLE)
  # nfd's macOS file dialog is Objective-C
  enable_language(OBJC)
endif()

# ---- Build type -----------------------------------------------------------------------------------
get_property(SCESS_MULTI_CONFIG GLOBAL PROPERTY GENERATOR_IS_MULTI_CONFIG)
set(SCESS_BUILD_TYPES Debug Release)
if(SCESS_MULTI_CONFIG)
  set(CMAKE_CONFIGURATION_TYPES ${SCESS_BUILD_TYPES} CACHE STRING "" FORCE)
else()
  if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release CACHE STRING "Build type" FORCE)
  endif()
  set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS ${SCESS_BUILD_TYPES})
  if(NOT CMAKE_BUILD_TYPE IN_LIST SCESS_BUILD_TYPES)
    message(FATAL_ERROR "CMAKE_BUILD_TYPE must be one of: ${SCESS_BUILD_TYPES} (got '${CMAKE_BUILD_TYPE}')")
  endif()
endif()

set(CMAKE_EXPORT_COMPILE_COMMANDS ON)
set(CMAKE_POSITION_INDEPENDENT_CODE ON)

set(SCESS_SRC_DIR ${PROJECT_SOURCE_DIR}/src)
set(SCESS_LIB_DIR ${PROJECT_SOURCE_DIR}/Libraries)
set(SCESS_ASSETS_DIR ${PROJECT_SOURCE_DIR}/Assets)

list(APPEND CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake)
include(SCESSHelpers)

# ---- Dependencies ---------------------------------------------------------------------------------
find_package(Boost CONFIG REQUIRED COMPONENTS filesystem)
find_package(fmt CONFIG REQUIRED)
find_package(Libint2 CONFIG REQUIRED)
find_package(Libxc CONFIG REQUIRED)
find_package(LAPACK REQUIRED)
find_library(LIB_LAPACKE lapacke REQUIRED)

# Dahlia is built as part of SCESS (in the same configuration), as its README recommends. This
# provides the Dahlia, DahliaTools and ImGui targets and dahlia_compile_shaders(). Dahlia needs
# glslc, GLFW, GLEW, GLM and FreeType (macOS: brew install shaderc glfw glew glm freetype).
get_filename_component(SCESS_DAHLIA_DEFAULT_DIR "${PROJECT_SOURCE_DIR}/../Dahlia" REALPATH)
set(SCESS_DAHLIA_DIR "${SCESS_DAHLIA_DEFAULT_DIR}" CACHE PATH "Dahlia graphics engine source checkout")
if(NOT EXISTS "${SCESS_DAHLIA_DIR}/src/Dahlia/CMakeLists.txt")
  message(FATAL_ERROR
    "\"${SCESS_DAHLIA_DIR}\" is not a Dahlia source checkout (src/Dahlia/CMakeLists.txt is missing). "
    "Clone Dahlia next to SCESS or specify one with -DSCESS_DAHLIA_DIR=...")
endif()
message(STATUS "Building Dahlia from: ${SCESS_DAHLIA_DIR}")
add_subdirectory("${SCESS_DAHLIA_DIR}" "${CMAKE_BINARY_DIR}/Dahlia")

# DahliaTools keeps its include directory private; QuantumSCESS and SCESS use its headers
# (Dahlia/Tools/*.h, Dahlia/Utilities/Memory.h) without linking libDahlia.
set(DAHLIA_INCLUDE_DIR "${SCESS_DAHLIA_DIR}/src")

# Git submodules: ImPlot, nfd and the header only libraries
add_subdirectory(Libraries)

# ---- Targets --------------------------------------------------------------------------------------
include(QuantumSCESS)
include(SCESS)
include(ViewerSCESS)
