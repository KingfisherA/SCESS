cmake_minimum_required(VERSION 3.16)

set(PREMAKE_BUILD_TYPES "Debug" "Release")
get_property(multi_config GLOBAL PROPERTY GENERATOR_IS_MULTI_CONFIG)
if(multi_config)
  set(CMAKE_CONFIGURATION_TYPES "${PREMAKE_BUILD_TYPES}" CACHE STRING "list of supported configuration types" FORCE)
else()
  set(CMAKE_BUILD_TYPE "Debug" CACHE STRING "Build Type of the project.")
  set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS "${PREMAKE_BUILD_TYPES}")
  if(NOT CMAKE_BUILD_TYPE IN_LIST PREMAKE_BUILD_TYPES)
    message(FATAL_ERROR
      "Invalid build type '${CMAKE_BUILD_TYPE}'.
      CMAKE_BUILD_TYPE must be any one of the possible values:
      ${PREMAKE_BUILD_TYPES}"
    )
  endif()
endif()

if(NOT IMGUI_INCLUDE_DIRS)
    message(FATAL_ERROR "ImGui include directory is not set. This is required for ImPlot to build.")
endif()

get_filename_component(IMPLOT_DIR "${CMAKE_CURRENT_SOURCE_DIR}/../ImPlot" ABSOLUTE)

project("ImPlot")
include(ImPlot.cmake)
