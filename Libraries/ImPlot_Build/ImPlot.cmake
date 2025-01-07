add_library("ImPlot" STATIC
    "../ImPlot/.github/example_implot.cpp"
   "../ImPlot/cmake-build-release/CMakeFiles/3.28.3/CompilerIdCXX/CMakeCXXCompilerId.cpp"
  "../ImPlot/implot.cpp"
  "../ImPlot/implot_items.cpp"
)
if(CMAKE_BUILD_TYPE STREQUAL Debug)
  set_target_properties("ImPlot" PROPERTIES
    OUTPUT_NAME "ImPlot"
    ARCHIVE_OUTPUT_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/../ImPlot/bin/linux/Debug-x64"
    LIBRARY_OUTPUT_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/../ImPlot/bin/linux/Debug-x64"
    RUNTIME_OUTPUT_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/../ImPlot/bin/linux/Debug-x64"
  )
endif()
target_include_directories("ImPlot" PRIVATE
  $<$<CONFIG:Debug>:${CMAKE_CURRENT_SOURCE_DIR}/../ImPlot>
  $<$<CONFIG:Debug>:${CMAKE_CURRENT_SOURCE_DIR}/../GraphicsEngine-1.0.0/include/ImGui>
)
target_compile_definitions("ImPlot" PRIVATE
  $<$<CONFIG:Debug>:CONFIG_DEBUG>
  $<$<CONFIG:Debug>:_DEBUG>
)
target_link_directories("ImPlot" PRIVATE
)
target_link_libraries("ImPlot"
)
target_compile_options("ImPlot" PRIVATE
  $<$<AND:$<CONFIG:Debug>,$<COMPILE_LANGUAGE:C>>:-m64>
  $<$<AND:$<CONFIG:Debug>,$<COMPILE_LANGUAGE:C>>:-O0>
  $<$<AND:$<CONFIG:Debug>,$<COMPILE_LANGUAGE:C>>:-g>
  $<$<AND:$<CONFIG:Debug>,$<COMPILE_LANGUAGE:CXX>>:-m64>
  $<$<AND:$<CONFIG:Debug>,$<COMPILE_LANGUAGE:CXX>>:-O0>
  $<$<AND:$<CONFIG:Debug>,$<COMPILE_LANGUAGE:CXX>>:-g>
  $<$<AND:$<CONFIG:Debug>,$<COMPILE_LANGUAGE:CXX>>:-std=c++20>
)
if(CMAKE_BUILD_TYPE STREQUAL Debug)
  set_target_properties("ImPlot" PROPERTIES
    CXX_STANDARD 20
    CXX_STANDARD_REQUIRED YES
    CXX_EXTENSIONS NO
    POSITION_INDEPENDENT_CODE False
    INTERPROCEDURAL_OPTIMIZATION False
  )
endif()
if(CMAKE_BUILD_TYPE STREQUAL Release)
  set_target_properties("ImPlot" PROPERTIES
    OUTPUT_NAME "ImPlot"
    ARCHIVE_OUTPUT_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/../ImPlot/bin/linux/Release-x64"
    LIBRARY_OUTPUT_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/../ImPlot/bin/linux/Release-x64"
    RUNTIME_OUTPUT_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/../ImPlot/bin/linux/Release-x64"
  )
endif()
target_include_directories("ImPlot" PRIVATE
  $<$<CONFIG:Release>:${CMAKE_CURRENT_SOURCE_DIR}/../ImPlot>
  $<$<CONFIG:Release>:${CMAKE_CURRENT_SOURCE_DIR}/../GraphicsEngine-1.0.0/include/ImGui>
)
target_compile_definitions("ImPlot" PRIVATE
  $<$<CONFIG:Release>:CONFIG_RELEASE>
)
target_link_directories("ImPlot" PRIVATE
)
target_link_libraries("ImPlot"
)
target_compile_options("ImPlot" PRIVATE
  $<$<AND:$<CONFIG:Release>,$<COMPILE_LANGUAGE:C>>:-m64>
  $<$<AND:$<CONFIG:Release>,$<COMPILE_LANGUAGE:C>>:-O3>
  $<$<AND:$<CONFIG:Release>,$<COMPILE_LANGUAGE:CXX>>:-m64>
  $<$<AND:$<CONFIG:Release>,$<COMPILE_LANGUAGE:CXX>>:-O3>
  $<$<AND:$<CONFIG:Release>,$<COMPILE_LANGUAGE:CXX>>:-std=c++20>
)
if(CMAKE_BUILD_TYPE STREQUAL Release)
  set_target_properties("ImPlot" PROPERTIES
    CXX_STANDARD 20
    CXX_STANDARD_REQUIRED YES
    CXX_EXTENSIONS NO
    POSITION_INDEPENDENT_CODE False
    INTERPROCEDURAL_OPTIMIZATION False
  )
endif()
