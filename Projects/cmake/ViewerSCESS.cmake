add_executable("ViewerSCESS"
    "src/Resources/AssetFolder.cpp"
    "src/Resources/Resources.h"
    "src/Resources/ShaderSources.cpp"
    "src/Utilities/Assets.h"
    "src/Utilities/Intrinsics.h"
    "src/Utilities/Utilities.cpp"
    "src/Utilities/Utilities.h"
      "src/Viewer/Application/Application.cpp"
      "src/Viewer/Application/Application.h"
      "src/Viewer/GUI/AtomDraw.cpp"
      "src/Viewer/GUI/AtomDraw.h"
      "src/Viewer/GUI/CalculationLayer.cpp"
      "src/Viewer/GUI/CalculationLayer.h"
      "src/Viewer/GUI/ImSpinner.h"
      "src/Viewer/Isosurfaces/IsosurfaceConstants.cpp"
      "src/Viewer/Isosurfaces/Isosurfaces.cpp"
      "src/Viewer/Isosurfaces/Isosurfaces.h"
)
if(CMAKE_BUILD_TYPE STREQUAL Debug)
  add_dependencies("ViewerSCESS"
    "QuantumSCESS"
  )
  set_target_properties("ViewerSCESS" PROPERTIES
    OUTPUT_NAME "ViewerSCESS"
    ARCHIVE_OUTPUT_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/bin/linux/Debug-x64"
    LIBRARY_OUTPUT_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/bin/linux/Debug-x64"
    RUNTIME_OUTPUT_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/bin/linux/Debug-x64"
  )
endif()
target_include_directories("ViewerSCESS" PRIVATE
  $<$<CONFIG:Debug>:${CMAKE_CURRENT_SOURCE_DIR}/src>
  $<$<CONFIG:Debug>:${CMAKE_CURRENT_SOURCE_DIR}/Libraries/GraphicsEngine-1.0.0/include>
  $<$<CONFIG:Debug>:${CMAKE_CURRENT_SOURCE_DIR}/Libraries/glm>
  $<$<CONFIG:Debug>:${CMAKE_CURRENT_SOURCE_DIR}/Libraries/cppitertools>
  $<$<CONFIG:Debug>:${CMAKE_CURRENT_SOURCE_DIR}/Libraries/JSON/include>
  $<$<CONFIG:Debug>:${CMAKE_CURRENT_SOURCE_DIR}/Libraries/Eigen>
  $<$<CONFIG:Debug>:${CMAKE_CURRENT_SOURCE_DIR}/Libraries/LibInt2/include>
  $<$<CONFIG:Debug>:${CMAKE_CURRENT_SOURCE_DIR}/Libraries/nfd/src/include>
  $<$<CONFIG:Debug>:${CMAKE_CURRENT_SOURCE_DIR}/Libraries/LibXC/lib/include>
)
target_compile_definitions("ViewerSCESS" PRIVATE
  $<$<CONFIG:Debug>:OS_LINUX>
  $<$<CONFIG:Debug>:CONFIG_DEBUG>
  $<$<CONFIG:Debug>:_DEBUG>
)
target_link_directories("ViewerSCESS" PRIVATE
  $<$<CONFIG:Debug>:${CMAKE_CURRENT_SOURCE_DIR}/Libraries/GraphicsEngine-1.0.0>
  $<$<CONFIG:Debug>:${CMAKE_CURRENT_SOURCE_DIR}/Libraries/nfd/build/lib/Release/x64>
  $<$<CONFIG:Debug>:${CMAKE_CURRENT_SOURCE_DIR}/Libraries/LibXC/lib>
)
target_link_libraries("ViewerSCESS"
  $<$<CONFIG:Debug>:QuantumSCESS>
  $<$<CONFIG:Debug>:GraphicsEngine>
  $<$<CONFIG:Debug>:xc.a>
  $<$<CONFIG:Debug>:nfd.a>
  $<$<CONFIG:Debug>:gtk-3>
  $<$<CONFIG:Debug>:gdk-3>
  $<$<CONFIG:Debug>:glib-2.0>
)
target_compile_options("ViewerSCESS" PRIVATE
  $<$<AND:$<CONFIG:Debug>,$<COMPILE_LANGUAGE:C>>:-m64>
  $<$<AND:$<CONFIG:Debug>,$<COMPILE_LANGUAGE:C>>:-O0>
  $<$<AND:$<CONFIG:Debug>,$<COMPILE_LANGUAGE:C>>:-g>
  $<$<AND:$<CONFIG:Debug>,$<COMPILE_LANGUAGE:CXX>>:-m64>
  $<$<AND:$<CONFIG:Debug>,$<COMPILE_LANGUAGE:CXX>>:-O0>
  $<$<AND:$<CONFIG:Debug>,$<COMPILE_LANGUAGE:CXX>>:-g>
  $<$<AND:$<CONFIG:Debug>,$<COMPILE_LANGUAGE:CXX>>:-std=c++20>
)
if(CMAKE_BUILD_TYPE STREQUAL Debug)
  set_target_properties("ViewerSCESS" PROPERTIES
    CXX_STANDARD 20
    CXX_STANDARD_REQUIRED YES
    CXX_EXTENSIONS NO
    POSITION_INDEPENDENT_CODE False
    INTERPROCEDURAL_OPTIMIZATION False
  )
endif()
if(CMAKE_BUILD_TYPE STREQUAL Debug)
target_precompile_headers("ViewerSCESS" PUBLIC ./src/pch.h)
endif()
if(CMAKE_BUILD_TYPE STREQUAL Release)
  add_dependencies("ViewerSCESS"
    "QuantumSCESS"
  )
  set_target_properties("ViewerSCESS" PROPERTIES
    OUTPUT_NAME "ViewerSCESS"
    ARCHIVE_OUTPUT_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/bin/linux/Release-x64"
    LIBRARY_OUTPUT_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/bin/linux/Release-x64"
    RUNTIME_OUTPUT_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/bin/linux/Release-x64"
  )
endif()
target_include_directories("ViewerSCESS" PRIVATE
  $<$<CONFIG:Release>:${CMAKE_CURRENT_SOURCE_DIR}/src>
  $<$<CONFIG:Release>:${CMAKE_CURRENT_SOURCE_DIR}/Libraries/GraphicsEngine-1.0.0/include>
  $<$<CONFIG:Release>:${CMAKE_CURRENT_SOURCE_DIR}/Libraries/glm>
  $<$<CONFIG:Release>:${CMAKE_CURRENT_SOURCE_DIR}/Libraries/cppitertools>
  $<$<CONFIG:Release>:${CMAKE_CURRENT_SOURCE_DIR}/Libraries/JSON/include>
  $<$<CONFIG:Release>:${CMAKE_CURRENT_SOURCE_DIR}/Libraries/Eigen>
  $<$<CONFIG:Release>:${CMAKE_CURRENT_SOURCE_DIR}/Libraries/LibInt2/include>
  $<$<CONFIG:Release>:${CMAKE_CURRENT_SOURCE_DIR}/Libraries/nfd/src/include>
  $<$<CONFIG:Release>:${CMAKE_CURRENT_SOURCE_DIR}/Libraries/LibXC/lib/include>
)
target_compile_definitions("ViewerSCESS" PRIVATE
  $<$<CONFIG:Release>:OS_LINUX>
  $<$<CONFIG:Release>:CONFIG_RELEASE>
)
target_link_directories("ViewerSCESS" PRIVATE
  $<$<CONFIG:Release>:${CMAKE_CURRENT_SOURCE_DIR}/Libraries/GraphicsEngine-1.0.0>
  $<$<CONFIG:Release>:${CMAKE_CURRENT_SOURCE_DIR}/Libraries/nfd/build/lib/Release/x64>
  $<$<CONFIG:Release>:${CMAKE_CURRENT_SOURCE_DIR}/Libraries/LibXC/lib>
)
target_link_libraries("ViewerSCESS"
  $<$<CONFIG:Release>:QuantumSCESS>
  $<$<CONFIG:Release>:GraphicsEngine>
  $<$<CONFIG:Release>:xc.a>
  $<$<CONFIG:Release>:nfd.a>
  $<$<CONFIG:Release>:gtk-3>
  $<$<CONFIG:Release>:gdk-3>
  $<$<CONFIG:Release>:glib-2.0>
)
target_compile_options("ViewerSCESS" PRIVATE
  $<$<AND:$<CONFIG:Release>,$<COMPILE_LANGUAGE:C>>:-m64>
  $<$<AND:$<CONFIG:Release>,$<COMPILE_LANGUAGE:C>>:-O3>
  $<$<AND:$<CONFIG:Release>,$<COMPILE_LANGUAGE:C>>:-mavx2>
  $<$<AND:$<CONFIG:Release>,$<COMPILE_LANGUAGE:CXX>>:-m64>
  $<$<AND:$<CONFIG:Release>,$<COMPILE_LANGUAGE:CXX>>:-O3>
  $<$<AND:$<CONFIG:Release>,$<COMPILE_LANGUAGE:CXX>>:-mavx2>
  $<$<AND:$<CONFIG:Release>,$<COMPILE_LANGUAGE:CXX>>:-std=c++20>
)
if(CMAKE_BUILD_TYPE STREQUAL Release)
  set_target_properties("ViewerSCESS" PROPERTIES
    CXX_STANDARD 20
    CXX_STANDARD_REQUIRED YES
    CXX_EXTENSIONS NO
    POSITION_INDEPENDENT_CODE False
    INTERPROCEDURAL_OPTIMIZATION False
  )
endif()
if(CMAKE_BUILD_TYPE STREQUAL Release)
target_precompile_headers("ViewerSCESS" PUBLIC ./src/pch.h)
endif()