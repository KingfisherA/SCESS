add_executable("ViewerSCESS"
      "src/Explorer/Application/Application.cpp"
      "src/Explorer/Application/Application.h"
      "src/Explorer/GUI/AtomStyle.cpp"
      "src/Explorer/GUI/AtomStyle.h"
      "src/Explorer/GUI/CalculationLayer.cpp"
      "src/Explorer/GUI/CalculationLayer.h"
      "src/Explorer/GUI/ImSpinner.h"
      "src/Explorer/GUI/ShaderSources.cpp"
      "src/Explorer/GUI/ShaderSources.h"
      "src/Explorer/Isosurfaces/IsosurfaceConstants.cpp"
      "src/Explorer/Isosurfaces/Isosurfaces.cpp"
      "src/Explorer/Isosurfaces/Isosurfaces.h"
    "src/Utilities/Intrinsics.h"
    "src/Utilities/Resources.h"
    "src/Utilities/Utilities.cpp"
    "src/Utilities/Utilities.h"
)
if(CMAKE_BUILD_TYPE STREQUAL Debug)
  add_dependencies("ViewerSCESS"
    "QuantumSCESS"
  )
  set_target_properties("ViewerSCESS" PROPERTIES
    OUTPUT_NAME "ViewerSCESS"
    ARCHIVE_OUTPUT_DIRECTORY "/media/charlie/Data/dev/C++/SCESS/bin/linux/Debug-x64"
    LIBRARY_OUTPUT_DIRECTORY "/media/charlie/Data/dev/C++/SCESS/bin/linux/Debug-x64"
    RUNTIME_OUTPUT_DIRECTORY "/media/charlie/Data/dev/C++/SCESS/bin/linux/Debug-x64"
  )
endif()
target_include_directories("ViewerSCESS" PRIVATE
  $<$<CONFIG:Debug>:/media/charlie/Data/dev/C++/SCESS/src>
  $<$<CONFIG:Debug>:/usr/include>
  $<$<CONFIG:Debug>:/usr/local/include>
  $<$<CONFIG:Debug>:/media/charlie/Data/dev/C++/GraphicsEngine/src>
  $<$<CONFIG:Debug>:/media/charlie/Data/dev/Libraries/json-develop/include>
  $<$<CONFIG:Debug>:/media/charlie/Data/dev/Libraries/eigen-3.4.90>
  $<$<CONFIG:Debug>:/media/charlie/Data/dev/Libraries/glm>
  $<$<CONFIG:Debug>:/media/charlie/Data/dev/Libraries>
  $<$<CONFIG:Debug>:/media/charlie/Data/dev/Libraries/ImGui>
  $<$<CONFIG:Debug>:/media/charlie/Data/dev/Libraries/nativefiledialog/src/include>
  $<$<CONFIG:Debug>:/opt/etsf/include>
)
target_compile_definitions("ViewerSCESS" PRIVATE
  $<$<CONFIG:Debug>:OS_LINUX>
  $<$<CONFIG:Debug>:CONFIG_DEBUG>
  $<$<CONFIG:Debug>:_DEBUG>
)
target_link_directories("ViewerSCESS" PRIVATE
  $<$<CONFIG:Debug>:/media/charlie/Data/dev/C++/GraphicsEngine/bin/linux/Release-x64/GraphicsEngine>
  $<$<CONFIG:Debug>:/media/charlie/Data/dev/Libraries/ImGui/bin/linux/Release-x64>
  $<$<CONFIG:Debug>:/media/charlie/Data/dev/Libraries/nativefiledialog/build/lib/Release/x64>
  $<$<CONFIG:Debug>:/opt/etsf/lib>
)
target_link_libraries("ViewerSCESS"
  $<$<CONFIG:Debug>:QuantumSCESS>
  $<$<CONFIG:Debug>:GraphicsEngine>
  $<$<CONFIG:Debug>:glfw3>
  $<$<CONFIG:Debug>:ImGui>
  $<$<CONFIG:Debug>:GLEW.a>
  $<$<CONFIG:Debug>:GL>
  $<$<CONFIG:Debug>:m>
  $<$<CONFIG:Debug>:Xrandr>
  $<$<CONFIG:Debug>:Xi>
  $<$<CONFIG:Debug>:X11>
  $<$<CONFIG:Debug>:Xxf86vm>
  $<$<CONFIG:Debug>:pthread>
  $<$<CONFIG:Debug>:nfd>
  $<$<CONFIG:Debug>:gtk-3>
  $<$<CONFIG:Debug>:z>
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
target_precompile_headers("ViewerSCESS" PUBLIC src/pch.h)
endif()
if(CMAKE_BUILD_TYPE STREQUAL Release)
  add_dependencies("ViewerSCESS"
    "QuantumSCESS"
  )
  set_target_properties("ViewerSCESS" PROPERTIES
    OUTPUT_NAME "ViewerSCESS"
    ARCHIVE_OUTPUT_DIRECTORY "/media/charlie/Data/dev/C++/SCESS/bin/linux/Release-x64"
    LIBRARY_OUTPUT_DIRECTORY "/media/charlie/Data/dev/C++/SCESS/bin/linux/Release-x64"
    RUNTIME_OUTPUT_DIRECTORY "/media/charlie/Data/dev/C++/SCESS/bin/linux/Release-x64"
  )
endif()
target_include_directories("ViewerSCESS" PRIVATE
  $<$<CONFIG:Release>:/media/charlie/Data/dev/C++/SCESS/src>
  $<$<CONFIG:Release>:/usr/include>
  $<$<CONFIG:Release>:/usr/local/include>
  $<$<CONFIG:Release>:/media/charlie/Data/dev/C++/GraphicsEngine/src>
  $<$<CONFIG:Release>:/media/charlie/Data/dev/Libraries/json-develop/include>
  $<$<CONFIG:Release>:/media/charlie/Data/dev/Libraries/eigen-3.4.90>
  $<$<CONFIG:Release>:/media/charlie/Data/dev/Libraries/glm>
  $<$<CONFIG:Release>:/media/charlie/Data/dev/Libraries>
  $<$<CONFIG:Release>:/media/charlie/Data/dev/Libraries/ImGui>
  $<$<CONFIG:Release>:/media/charlie/Data/dev/Libraries/nativefiledialog/src/include>
  $<$<CONFIG:Release>:/opt/etsf/include>
)
target_compile_definitions("ViewerSCESS" PRIVATE
  $<$<CONFIG:Release>:OS_LINUX>
  $<$<CONFIG:Release>:CONFIG_RELEASE>
)
target_link_directories("ViewerSCESS" PRIVATE
  $<$<CONFIG:Release>:/media/charlie/Data/dev/C++/GraphicsEngine/bin/linux/Release-x64/GraphicsEngine>
  $<$<CONFIG:Release>:/media/charlie/Data/dev/Libraries/ImGui/bin/linux/Release-x64>
  $<$<CONFIG:Release>:/media/charlie/Data/dev/Libraries/nativefiledialog/build/lib/Release/x64>
  $<$<CONFIG:Release>:/opt/etsf/lib>
)
target_link_libraries("ViewerSCESS"
  $<$<CONFIG:Release>:QuantumSCESS>
  $<$<CONFIG:Release>:GraphicsEngine>
  $<$<CONFIG:Release>:glfw3>
  $<$<CONFIG:Release>:ImGui>
  $<$<CONFIG:Release>:GLEW.a>
  $<$<CONFIG:Release>:GL>
  $<$<CONFIG:Release>:m>
  $<$<CONFIG:Release>:Xrandr>
  $<$<CONFIG:Release>:Xi>
  $<$<CONFIG:Release>:X11>
  $<$<CONFIG:Release>:Xxf86vm>
  $<$<CONFIG:Release>:pthread>
  $<$<CONFIG:Release>:nfd>
  $<$<CONFIG:Release>:gtk-3>
  $<$<CONFIG:Release>:z>
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
target_precompile_headers("ViewerSCESS" PUBLIC src/pch.h)
endif()