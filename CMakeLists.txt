# CMake Build Script for the SCESS project
# Includes the projectes for the quantum shared libraries, the GUI and the command line executable of SCESS 

cmake_minimum_required(VERSION 3.16)

set(ALLOWED_BUILD_TYPES "Debug" "Release")
if(NOT CMAKE_BUILD_TYPE IN_LIST ALLOWED_BUILD_TYPES)
  message(FATAL_ERROR
    "Invalid build type '${CMAKE_BUILD_TYPE}'.
    CMAKE_BUILD_TYPE must be any one of the possible values:
    ${ALLOWED_BUILD_TYPES}"
  )
endif()

if(NOT CMAKE_SYSTEM_NAME)
  set(CMAKE_SYSTEM_NAME "Linux")
endif()

if(NOT SCESS_GENG_DIR)
    get_filename_component(SCESS_GENG_DIR "../GraphicsEngine" ABSOLUTE)
    message(NOTICE "Graphics Engine directory not set. Setting to default: \"../GraphicsEngine\"")
endif()

get_filename_component(SCESS_DIR ${CMAKE_CURRENT_SOURCE_DIR} ABSOLUTE)

set(SCESS_PLATFORM ${CMAKE_SYSTEM_NAME})
string(TOLOWER ${SCESS_PLATFORM} SCESS_PLATFORM)
set(SCESS_SRC_DIR ${SCESS_DIR}/src)
set(SCESS_LIB_DIR ${SCESS_DIR}/Libraries)
set(SCESS_BIN_DIR ${SCESS_DIR}/bin/${SCESS_PLATFORM})

project("SCESS" LANGUAGES CXX)

find_package(Boost CONFIG REQUIRED COMPONENTS filesystem)
find_package(Libxc REQUIRED)
find_package(libint2 REQUIRED)
find_library(GENG_LIBRARIES GraphicsEngine REQUIRED
  PATHS "${SCESS_GENG_DIR}/bin"
  PATH_SUFFIXES "linux/Release-x64")

if(NOT GENG_LIBRARIES STREQUAL "GENG_LIBRARIES-NOTFOUND")
    set(IMGUI_INCLUDE_DIRS ${SCESS_GENG_DIR}/include/ImGui)
    set(GENG_INCLUDE_DIRS ${SCESS_GENG_DIR}/include)
endif()

include(Projects/cmake/QuantumSCESS.cmake)
include(Projects/cmake/SCESS.cmake)
include(Projects/cmake/ViewerSCESS.cmake)
