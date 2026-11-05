# Third party libraries vendored as git submodules. Eigen, cppitertools and JSON are header only
# and are added to the include paths of the SCESS targets; ImPlot and nfd are built here.

foreach(FILE Eigen/Eigen cppitertools/cppitertools/range.hpp JSON/include/nlohmann/json.hpp ImPlot/implot.h nfd/src/include/nfd.h)
  if(NOT EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/${FILE})
    message(FATAL_ERROR
      "Libraries/${FILE} is missing. Run: git submodule update --init --recursive")
  endif()
endforeach()

# ---- ImPlot ---------------------------------------------------------------------------------------
# Only ImGui's headers are needed: its symbols are resolved from libDahlia at link time, so ImPlot
# shares Dahlia's ImGui context rather than getting a second copy of ImGui.
set(IMPLOT_DIR ${CMAKE_CURRENT_SOURCE_DIR}/ImPlot)
add_library(ImPlot STATIC
  ${IMPLOT_DIR}/implot.cpp
  ${IMPLOT_DIR}/implot_items.cpp
)
target_include_directories(ImPlot
  PUBLIC ${IMPLOT_DIR}
  PRIVATE $<TARGET_PROPERTY:ImGui,INTERFACE_INCLUDE_DIRECTORIES>
)
set_target_properties(ImPlot PROPERTIES
  CXX_STANDARD 20
  CXX_STANDARD_REQUIRED ON
  CXX_EXTENSIONS OFF
)

# ---- nfd (Native File Dialog) ---------------------------------------------------------------------
set(NFD_DIR ${CMAKE_CURRENT_SOURCE_DIR}/nfd/src)
add_library(nfd STATIC ${NFD_DIR}/nfd_common.c)
target_include_directories(nfd PUBLIC ${NFD_DIR}/include)

if(APPLE)
  target_sources(nfd PRIVATE ${NFD_DIR}/nfd_cocoa.m)
  # Upstream nfd (unmaintained) still uses -[NSSavePanel setAllowedFileTypes:], deprecated in macOS 12
  set_source_files_properties(${NFD_DIR}/nfd_cocoa.m PROPERTIES COMPILE_OPTIONS -Wno-deprecated-declarations)
  target_link_libraries(nfd PUBLIC $<LINK_LIBRARY:FRAMEWORK,AppKit>)
elseif(WIN32)
  target_sources(nfd PRIVATE ${NFD_DIR}/nfd_win.cpp)
  target_link_libraries(nfd PUBLIC comctl32 ole32 uuid)
else()
  find_package(PkgConfig REQUIRED)
  pkg_check_modules(GTK3 REQUIRED IMPORTED_TARGET gtk+-3.0)
  target_sources(nfd PRIVATE ${NFD_DIR}/nfd_gtk.c)
  target_link_libraries(nfd PUBLIC PkgConfig::GTK3)
endif()
