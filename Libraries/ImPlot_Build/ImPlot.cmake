add_library("ImPlot" STATIC
  "${IMPLOT_DIR}/implot.cpp"
  "${IMPLOT_DIR}/implot_items.cpp"
)

target_include_directories("ImPlot" PRIVATE
  ${IMPLOT_DIR}
  ${IMGUI_INCLUDE_DIRS}
)

target_compile_definitions("ImPlot" PRIVATE
  $<$<CONFIG:Release>:CONFIG_RELEASE>
  $<$<CONFIG:Debug>:CONFIG_DEBUG>
  $<$<CONFIG:Debug>:_DEBUG>
)

target_compile_options("ImPlot" PRIVATE
  -std=c++20
)

set_target_properties("ImPlot" PROPERTIES
  CXX_STANDARD 20
  CXX_STANDARD_REQUIRED YES
  CXX_EXTENSIONS NO
  POSITION_INDEPENDENT_CODE False
  INTERPROCEDURAL_OPTIMIZATION False
  OUTPUT_NAME "ImPlot"
)

if(CMAKE_BUILD_TYPE STREQUAL Debug)
  set_target_properties("ImPlot" PROPERTIES
      ARCHIVE_OUTPUT_DIRECTORY "${IMPLOT_DIR}/bin/Debug-x64"
      LIBRARY_OUTPUT_DIRECTORY "${IMPLOT_DIR}/bin/Debug-x64"
      RUNTIME_OUTPUT_DIRECTORY "${IMPLOT_DIR}/bin/Debug-x64"
  )

  target_compile_options("ImPlot" PRIVATE
    -O0
    -g
  )
endif()

if(CMAKE_BUILD_TYPE STREQUAL Release)
  set_target_properties("ImPlot" PROPERTIES
      ARCHIVE_OUTPUT_DIRECTORY "${IMPLOT_DIR}/bin/Release-x64"
      LIBRARY_OUTPUT_DIRECTORY "${IMPLOT_DIR}/bin/Release-x64"
      RUNTIME_OUTPUT_DIRECTORY "${IMPLOT_DIR}/bin/Release-x64"
  )

  target_compile_options("ImPlot" PRIVATE
    -O3
  )
endif()
