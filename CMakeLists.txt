# CMake Build Script for the SCESS project
# Includes the projectes for the quantum shared libraries, the GUI and the command line executable of SCESS 

cmake_minimum_required(VERSION 3.16)

set(SCESS_BUILD_TYPES "Debug" "Release")
if(CMAKE_CONFIGURATION_TYPES)
#  set(CMAKE_CONFIGURATION_TYPES SCESS_BUILD_TYPES)
  message(NOTICE "---Multiconfiguration build generator selected---")
elseif(NOT CMAKE_BUILD_TYPE IN_LIST SCESS_BUILD_TYPES)
  message(FATAL_ERROR
    "Invalid build type '${CMAKE_BUILD_TYPE}'.
    CMAKE_BUILD_TYPE must be any one of the possible values:
    ${SCESS_BUILD_TYPES}"
  )
endif()

project("SCESS" LANGUAGES CXX)

if(NOT SCESS_GENG_DIR)
    message(NOTICE "-- Graphics Engine directory not set. Setting to default: \"../GraphicsEngine\"")
    get_filename_component(SCESS_GENG_DIR "../GraphicsEngine" REALPATH)
    if(NOT (EXISTS "${SCESS_GENG_DIR}"))   
        message(FATAL_ERROR "!! Cannot find GraphicsEngine directory. Please specify with -DSCESS_GENG_DIR=...")
    else()
        message(NOTICE "-- GraphicsEngine directory successfully bound!")
    endif()
endif()

get_filename_component(SCESS_DIR ${CMAKE_CURRENT_SOURCE_DIR} ABSOLUTE)

set(SCESS_PLATFORM ${CMAKE_SYSTEM_NAME})
string(TOLOWER ${SCESS_PLATFORM} SCESS_PLATFORM)
set(SCESS_SRC_DIR ${SCESS_DIR}/src)
set(SCESS_LIB_DIR ${SCESS_DIR}/Libraries)
set(SCESS_BIN_DIR ${SCESS_DIR}/bin/${SCESS_PLATFORM})

find_package(Boost CONFIG COMPONENTS filesystem locale math preprocessor)

if(NOT Boost_FOUND)
    find_package(Boost COMPONENTS REQUIRED filesystem locale)
endif()

find_package(Libxc CONFIG REQUIRED)
find_package(libint2 CONFIG REQUIRED)

find_library(GENG_LIBRARIES GraphicsEngine REQUIRED
  PATHS "${SCESS_GENG_DIR}/bin"
  PATH_SUFFIXES "${SCESS_PLATFORM}/Release")
if(NOT (GENG_LIBRARIES STREQUAL "GENG_LIBRARIES-NOTFOUND"))
    get_filename_component(IMGUI_INCLUDE_DIRS "${SCESS_GENG_DIR}/include/ImGui" REALPATH)
    get_filename_component(GENG_INCLUDE_DIRS "${SCESS_GENG_DIR}/include" REALPATH)

    if((NOT EXISTS "${IMGUI_INCLUDE_DIRS}") OR (NOT EXISTS "${GENG_INCLUDE_DIRS}"))   
        message(FATAL_ERROR "!! Cannot find GraphicsEngine include directory and or the GraphicsEngine's ImGui include directory. Please ensure that the include files for ImGui and GraphicsEngine are placed in a directory named 'include' in the folder:  \"${SCESS_GENG_DIR}\" ")
    else()
        message(NOTICE "-- Found GraphicsEngine include directories: ${GENG_INCLUDE_DIRS}, ${IMGUI_INCLUDE_DIRS}")
    endif()

    message("-- Found GraphicsEngine: ${GENG_LIBRARIES}")
endif()

find_library(GENG_TOOLS_LIBRARIES GraphicsEngineTools REQUIRED
        PATHS "${SCESS_GENG_DIR}/bin"
        PATH_SUFFIXES "${SCESS_PLATFORM}/Release")
if(NOT (GENG_TOOLS_LIBRARIES STREQUAL "GENG_TOOLS_LIBRARIES-NOTFOUND"))
    if(NOT DEFINED GENG_INCLUDE_DIRS)
        get_filename_component(GENG_INCLUDE_DIRS "${SCESS_GENG_DIR}/include" REALPATH)
    endif()

    if((NOT EXISTS "${GENG_INCLUDE_DIRS}"))
        message(FATAL_ERROR "!! Cannot find GraphicsEngineTool's include directory !!")
    endif()

    message("-- Found GraphicsEngineTools: ${GENG_TOOLS_LIBRARIES}")
endif()

include(Projects/cmake/QuantumSCESS.cmake)
include(Projects/cmake/SCESS.cmake)
include(Projects/cmake/ViewerSCESS.cmake)
