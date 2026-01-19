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

if(NOT SCESS_DAHLIA_DIR)
    message(NOTICE "-- Dahlia Engine directory not set. Setting to default: \"../Dahlia\"")
    get_filename_component(SCESS_DAHLIA_DIR "../Dahlia" REALPATH)
    if(NOT (EXISTS "${SCESS_DAHLIA_DIR}"))   
        message(FATAL_ERROR "!! Cannot find Dahlia directory. Please specify with -DSCESS_DAHLIA_DIR=...")
    else()
        message(NOTICE "-- Dahlia Engine directory successfully bound!")
    endif()
endif()

get_filename_component(SCESS_DIR ${CMAKE_CURRENT_SOURCE_DIR} ABSOLUTE)

set(SCESS_PLATFORM ${CMAKE_SYSTEM_NAME})
string(TOLOWER ${SCESS_PLATFORM} SCESS_PLATFORM)
set(SCESS_SRC_DIR ${SCESS_DIR}/src)
set(SCESS_LIB_DIR ${SCESS_DIR}/Libraries)
set(SCESS_BIN_DIR ${SCESS_DIR}/bin/${SCESS_PLATFORM})

find_package(Boost CONFIG COMPONENTS filesystem locale math preprocessor)
find_package(Libint2 CONFIG REQUIRED)

if(NOT Boost_FOUND)
    find_package(Boost COMPONENTS REQUIRED filesystem locale)
endif()

find_package(fmt CONFIG REQUIRED)
find_package(Libxc CONFIG REQUIRED)

find_library(DAHLIA_LIBRARIES Dahlia REQUIRED
  PATHS "${SCESS_DAHLIA_DIR}/bin"
  PATH_SUFFIXES "${SCESS_PLATFORM}/Release")
if(NOT (DAHLIA_LIBRARIES STREQUAL "DAHLIA_LIBRARIES-NOTFOUND"))
    get_filename_component(IMGUI_INCLUDE_DIRS "${SCESS_DAHLIA_DIR}/include/ImGui" REALPATH)
    get_filename_component(DAHLIA_INCLUDE_DIRS "${SCESS_DAHLIA_DIR}/include" REALPATH)

    if((NOT EXISTS "${IMGUI_INCLUDE_DIRS}") OR (NOT EXISTS "${DAHLIA_INCLUDE_DIRS}"))   
        message(FATAL_ERROR "!! Cannot find Dahlia Engine include directory and or the Dahlie Engine's ImGui include directory. Please ensure that the include files for ImGui and Dahlia Engine are placed in a directory named 'include' in the folder:  \"${SCESS_DAHLIA_DIR}\" ")
    else()
        message(NOTICE "-- Found Dahlia Engine include directories: ${DAHLIA_INCLUDE_DIRS}, ${IMGUI_INCLUDE_DIRS}")
    endif()

    message("-- Found Dahlia Engine: ${DAHLIA_LIBRARIES}")
endif()

find_library(DHALIA_TOOLS_LIBRARIES DahliaTools REQUIRED
        PATHS "${SCESS_DAHLIA_DIR}/bin"
        PATH_SUFFIXES "${SCESS_PLATFORM}/Release")
if(NOT (DHALIA_TOOLS_LIBRARIES STREQUAL "DHALIA_TOOLS_LIBRARIES-NOTFOUND"))
    if(NOT DEFINED DAHLIA_INCLUDE_DIRS)
        get_filename_component(DAHLIA_INCLUDE_DIRS "${SCESS_DAHLIA_DIR}/include" REALPATH)
    endif()

    if((NOT EXISTS "${DAHLIA_INCLUDE_DIRS}"))
        message(FATAL_ERROR "!! Cannot find Dahlia Engine Tool's include directory !!")
    endif()

    message("-- Found Dahlia EngineTools: ${DHALIA_TOOLS_LIBRARIES}")
endif()

include(Projects/cmake/QuantumSCESS.cmake)
include(Projects/cmake/SCESS.cmake)
include(Projects/cmake/ViewerSCESS.cmake)
