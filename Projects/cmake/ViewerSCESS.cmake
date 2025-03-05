#Project for defining the OpenGL GUI app for SCESS
add_executable("ViewerSCESS"
  "${SCESS_SRC_DIR}/Resources/AssetFolder.cpp"
  "${SCESS_SRC_DIR}/Resources/AssetFolder.h"
  "${SCESS_SRC_DIR}/Resources/ShaderSources.cpp"
  "${SCESS_SRC_DIR}/Resources/ShaderSources.h"
  "${SCESS_SRC_DIR}/Utilities/Assets.h"
  "${SCESS_SRC_DIR}/Utilities/Intrinsics.h"
  "${SCESS_SRC_DIR}/Utilities/Utilities.cpp"
  "${SCESS_SRC_DIR}/Utilities/Utilities.h"
  "${SCESS_SRC_DIR}/Viewer/Application/Application.cpp"
  "${SCESS_SRC_DIR}/Viewer/Application/Application.h"
  "${SCESS_SRC_DIR}/Viewer/GUI/AtomDraw.cpp"
  "${SCESS_SRC_DIR}/Viewer/GUI/AtomDraw.h"
  "${SCESS_SRC_DIR}/Viewer/GUI/CalculationLayer.cpp"
  "${SCESS_SRC_DIR}/Viewer/GUI/CalculationLayer.h"
  "${SCESS_SRC_DIR}/Viewer/GUI/ImSpinner.h"
  "${SCESS_SRC_DIR}/Viewer/Isosurfaces/IsosurfaceConstants.cpp"
  "${SCESS_SRC_DIR}/Viewer/Isosurfaces/Isosurfaces.cpp"
  "${SCESS_SRC_DIR}/Viewer/Isosurfaces/Isosurfaces.h"
)

add_subdirectory(${SCESS_LIB_DIR}/ImPlot_Build)

add_dependencies("ViewerSCESS"
   "QuantumSCESS"
   "ImPlot"
)

if(CMAKE_SYSTEM_NAME STREQUAL Linux)

  find_package (PkgConfig REQUIRED)
  pkg_check_modules (GTK3 REQUIRED gtk+-3.0)
  target_include_directories("ViewerSCESS" PRIVATE ${GTK3_INCLUDE_DIRS})
  target_link_directories("ViewerSCESS" PRIVATE ${GTK3_LIBRARY_DIRS})
  target_link_libraries("ViewerSCESS" PRIVATE ${GTK3_LIBRARIES})

elseif(CMAKE_SYSTEM_NAME STREQUAL Windows)
  
  target_link_libraries("ViewerSCESS" PRIVATE comctl32.lib)

endif()


find_library(LIB_nfd nfd REQUIRED
        PATHS "${SCESS_LIB_DIR}/nfd/build"
        PATH_SUFFIXES "lib/Release/x64")

target_include_directories("ViewerSCESS" PRIVATE
  ${SCESS_SRC_DIR}
  ${SCESS_LIB_DIR}/ImPlot
  ${GENG_INCLUDE_DIRS}
  ${IMGUI_INCLUDE_DIRS}
  ${SCESS_LIB_DIR}/cppitertools
  ${SCESS_LIB_DIR}/JSON/include
  ${SCESS_LIB_DIR}/Eigen
  ${SCESS_LIB_DIR}/nfd/src/include
)

target_link_libraries("ViewerSCESS" PRIVATE
  fmt
  QuantumSCESS
  ${GENG_LIBRARIES}
  ImPlot
  ${LIB_nfd}
  Boost::filesystem
  $<TARGET_NAME_IF_EXISTS:Libint2::int2>
  $<TARGET_NAME_IF_EXISTS:Libint2::libint2>
  Libxc::xc
)

target_compile_definitions("ViewerSCESS" PRIVATE
        LIBINT2_CPLUSPLUS_STD=2020
        $<$<PLATFORM_ID:Linux>:OS_LINUX>
        $<$<PLATFORM_ID:Windows>:OS_WINDOWS>
        $<$<CONFIG:DEBUG>:CONFIG_DEBUG>
        $<$<CONFIG:DEBUG>:_DEBUG>
        $<$<CONFIG:RELEASE>:CONFIG_RELEASE>)

#target_precompile_headers("ViewerSCESS" PUBLIC ./src/Q_PCH.h)

set(COMPILE_OPTIONS)
if((CMAKE_CXX_COMPILER_ID STREQUAL Clang) OR 
   (CMAKE_CXX_COMPILER_ID STREQUAL GNU))
    list(APPEND COMPILE_OPTIONS
      -m64
      -fPIC
      -std=c++20
    )
    if(CMAKE_BUILD_TYPE STREQUAL Debug)
      list(APPEND COMPILE_OPTIONS
        -g
        -O0
      )
    elseif(CMAKE_BUILD_TYPE STREQUAL Release)
      list(APPEND COMPILE_OPTIONS
        -O3
        -mavx2
      )
    endif()
elseif(CMAKE_CXX_COMPILER_ID STREQUAL MSVC)
  list(APPEND COMPILE_OPTIONS
      /std:c++20
    )
    if(CMAKE_BUILD_TYPE STREQUAL Debug)
      list(APPEND COMPILE_OPTIONS
        /Zi
        /Od
      )
    elseif(CMAKE_BUILD_TYPE STREQUAL Release)
      list(APPEND COMPILE_OPTIONS /O2)
    endif()
endif()
target_compile_options("ViewerSCESS" PRIVATE ${COMPILE_OPTIONS})

set_target_properties("ViewerSCESS" PROPERTIES
    CXX_STANDARD 20
    CXX_STANDARD_REQUIRED YES
    CXX_EXTENSIONS NO
    POSITION_INDEPENDENT_CODE True
    INTERPROCEDURAL_OPTIMIZATION False
    MSVC_RUNTIME_LIBRARY MultiThreaded$<$<CONFIG:Debug>:Debug>DLL
    OUTPUT_NAME "ViewerSCESS")

if(CMAKE_CONFIGURATION_TYPES)
  set(DEBUG_BIN_DIR "${SCESS_BIN_DIR}")
  set(RELEASE_BIN_DIR "${SCESS_BIN_DIR}")
else()
  set(DEBUG_BIN_DIR "${SCESS_BIN_DIR}/Debug")
  set(RELEASE_BIN_DIR "${SCESS_BIN_DIR}/Release")
endif()

if(CMAKE_BUILD_TYPE STREQUAL Debug)
  set_target_properties("ViewerSCESS" PROPERTIES
    ARCHIVE_OUTPUT_DIRECTORY "${DEBUG_BIN_DIR}"
    LIBRARY_OUTPUT_DIRECTORY "${DEBUG_BIN_DIR}"
    RUNTIME_OUTPUT_DIRECTORY "${DEBUG_BIN_DIR}"
  )
elseif(CMAKE_BUILD_TYPE STREQUAL Release)
  set_target_properties("ViewerSCESS" PROPERTIES
    ARCHIVE_OUTPUT_DIRECTORY "${RELEASE_BIN_DIR}"
    LIBRARY_OUTPUT_DIRECTORY "${RELEASE_BIN_DIR}"
    RUNTIME_OUTPUT_DIRECTORY "${RELEASE_BIN_DIR}"
  )
endif()
