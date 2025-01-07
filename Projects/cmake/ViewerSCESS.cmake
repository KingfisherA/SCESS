#Project for defining the OpenGL GUI app for SCESS
add_executable("ViewerSCESS"
  "${SCESS_SRC_DIR}/Resources/AssetFolder.cpp"
  "${SCESS_SRC_DIR}/Resources/Resources.h"
  "${SCESS_SRC_DIR}/Resources/ShaderSources.cpp"
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
target_compile_options("ViewerSCESS" PRIVATE
  -m64
  -std=c++20
)

find_package (PkgConfig REQUIRED)
pkg_check_modules (GTK3 REQUIRED gtk+-3.0)
find_library(LIB_nfd nfd REQUIRED
        PATHS "${SCESS_LIB_DIR}/nfd/build"
        PATH_SUFFIXES "lib/Release/x64")

target_include_directories("ViewerSCESS" PRIVATE
  ${SCESS_SRC_DIR}
  ${SCESS_LIB_DIR}/ImPlot
  ${SCESS_LIB_DIR}/GraphicsEngine-1.0.0/include
  ${SCESS_LIB_DIR}/GraphicsEngine-1.0.0/include/ImGui
  ${SCESS_LIB_DIR}/cppitertools
  ${SCESS_LIB_DIR}/JSON/include
  ${SCESS_LIB_DIR}/Eigen
  ${SCESS_LIB_DIR}/nfd/src/include
  ${GTK3_INCLUDE_DIRS}
)
target_link_directories("ViewerSCESS" PRIVATE
        ${GTK3_LIBRARY_DIRS})
target_link_libraries("ViewerSCESS" PRIVATE
  QuantumSCESS
  ${LIB_GraphicsEngine}
  ImPlot
  ${LIB_nfd}
  Boost::filesystem
  Libint2::libint2
  Libxc::xc
  ${GTK3_LIBRARIES}
)

target_compile_definitions("ViewerSCESS" PRIVATE
        OS_LINUX
        $<$<CONFIG:DEBUG>:CONFIG_DEBUG>
        $<$<CONFIG:DEBUG>:_DEBUG>
        $<$<CONFIG:RELEASE>:CONFIG_RELEASE>
)
target_precompile_headers("ViewerSCESS" PUBLIC ./src/pch.h)

if(CMAKE_BUILD_TYPE STREQUAL Debug)
  set_target_properties("ViewerSCESS" PROPERTIES
    CXX_STANDARD 20
    CXX_STANDARD_REQUIRED YES
    CXX_EXTENSIONS NO
    POSITION_INDEPENDENT_CODE False
    INTERPROCEDURAL_OPTIMIZATION False
    INTERFACE_LINK_LIBRARIES ""
    OUTPUT_NAME "ViewerSCESS"
    ARCHIVE_OUTPUT_DIRECTORY "${SCESS_BIN_DIR}/Debug-x64"
    LIBRARY_OUTPUT_DIRECTORY "${SCESS_BIN_DIR}/Debug-x64"
    RUNTIME_OUTPUT_DIRECTORY "${SCESS_BIN_DIR}/Debug-x64"
  )
  target_compile_options("ViewerSCESS" PRIVATE
    -O0
    -g
  ) 
endif()

if(CMAKE_BUILD_TYPE STREQUAL Release)
  set_target_properties("ViewerSCESS" PROPERTIES
    CXX_STANDARD 20
    CXX_STANDARD_REQUIRED YES
    CXX_EXTENSIONS NO
    POSITION_INDEPENDENT_CODE False
    INTERPROCEDURAL_OPTIMIZATION False
    OUTPUT_NAME "ViewerSCESS"
    ARCHIVE_OUTPUT_DIRECTORY "${SCESS_BIN_DIR}/Release-x64"
    LIBRARY_OUTPUT_DIRECTORY "${SCESS_BIN_DIR}/Release-x64"
    RUNTIME_OUTPUT_DIRECTORY "${SCESS_BIN_DIR}/Release-x64"
  )
  target_compile_options("ViewerSCESS" PRIVATE
    -O2
    -mavx2
  ) 
endif()