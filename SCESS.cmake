add_executable("SCESS"
  "SCESS.APS"
  "SCESS.rc"
  "resource.h"
      "src/3D/Buffers.cpp"
      "src/3D/Buffers.h"
      "src/3D/Framebuffer.cpp"
      "src/3D/Framebuffer.h"
      "src/3D/ShaderProgram.cpp"
      "src/3D/ShaderProgram.h"
      "src/3D/Texture.cpp"
      "src/3D/Texture.h"
      "src/3D/VertexArray.cpp"
      "src/3D/VertexArray.h"
        "src/3D/stb/stb_image.cpp"
        "src/3D/stb/stb_image.h"
      "src/Application/Application.cpp"
      "src/Application/Application.h"
      "src/Chemistry/PeriodicTable.cpp"
      "src/Chemistry/PeriodicTable.h"
    "src/Display.cpp"
    "src/Display.h"
      "src/Events/AppEvents.h"
      "src/Events/Events.h"
      "src/Events/KeyEvents.h"
      "src/Events/MouseEvents.h"
      "src/Events/WindowEvents.h"
      "src/GUI/CalculationLayer.cpp"
      "src/GUI/CalculationLayer.h"
      "src/GUI/ImGuiBuild.cpp"
      "src/GUI/ImGuiLayer.cpp"
      "src/GUI/ImGuiLayer.h"
      "src/GUI/ImSpinner.h"
      "src/Layers/Layer.cpp"
      "src/Layers/Layer.h"
      "src/Layers/LayerStack.cpp"
      "src/Layers/LayerStack.h"
      "src/Maths/Geometry.cpp"
      "src/Maths/Geometry.hpp"
      "src/Maths/Isosurfaces.cpp"
      "src/Maths/Isosurfaces.h"
      "src/Maths/Maths.cpp"
      "src/Maths/Maths.hpp"
      "src/Maths/Quadrature.cpp"
      "src/Maths/Quadrature.h"
      "src/Platform/API Core.h"
        "src/Platform/GLFW/GLFWDisplay.cpp"
        "src/Platform/GLFW/GLFWDisplay.h"
      "src/Platform/Intrinsics.h"
        "src/Platform/OpenGL/OpenGLBuffers.cpp"
        "src/Platform/OpenGL/OpenGLBuffers.h"
        "src/Platform/OpenGL/OpenGLCore.cpp"
        "src/Platform/OpenGL/OpenGLCore.h"
        "src/Platform/OpenGL/OpenGLFramebuffer.cpp"
        "src/Platform/OpenGL/OpenGLFramebuffer.h"
        "src/Platform/OpenGL/OpenGLShaderProgram.cpp"
        "src/Platform/OpenGL/OpenGLShaderProgram.h"
        "src/Platform/OpenGL/OpenGLTexture.cpp"
        "src/Platform/OpenGL/OpenGLTexture.h"
        "src/Platform/OpenGL/OpenGLVertexArray.cpp"
        "src/Platform/OpenGL/OpenGLVertexArray.h"
        "src/Platform/Windows/Console.cpp"
        "src/Platform/Windows/Console.h"
        "src/Quantum/Basis/BasisUtils.cpp"
        "src/Quantum/Basis/BasisUtils.h"
        "src/Quantum/Basis/Gaussian.cpp"
        "src/Quantum/Basis/Gaussian.h"
        "src/Quantum/DFT/DFT.cpp"
        "src/Quantum/DFT/DFT.h"
        "src/Quantum/DFT/DFTGrid.cpp"
        "src/Quantum/DFT/DFTGrid.h"
        "src/Quantum/DFT/XCFunctional.cpp"
        "src/Quantum/DFT/XCFunctional.h"
        "src/Quantum/HF/HartreeFock.cpp"
        "src/Quantum/HF/HartreeFock.h"
        "src/Quantum/Integrals/Gamma.cpp"
        "src/Quantum/Integrals/Gamma.h"
        "src/Quantum/Integrals/IntegralHelper.cpp"
        "src/Quantum/Integrals/IntegralHelper.h"
      "src/Quantum/Molecule.cpp"
      "src/Quantum/Molecule.h"
        "src/Quantum/Tests/Tests.h"
      "src/Scene/Camera.cpp"
      "src/Scene/Camera.hpp"
      "src/Utilities/Assert.h"
      "src/Utilities/IO.cpp"
      "src/Utilities/IO.h"
      "src/Utilities/Memory.h"
      "src/Utilities/PathUtils.h"
      "src/Utilities/StringOperations.cpp"
      "src/Utilities/StringOperations.h"
      "src/Utilities/Threading.h"
    "src/pch.cpp"
    "src/pch.h"
)
if(CMAKE_BUILD_TYPE STREQUAL Debug)
  set_target_properties("SCESS" PROPERTIES
    OUTPUT_NAME "SCESS"
    ARCHIVE_OUTPUT_DIRECTORY "/media/charlie/Data/dev/C++/SCESS/bin/linux/Debug-x64"
    LIBRARY_OUTPUT_DIRECTORY "/media/charlie/Data/dev/C++/SCESS/bin/linux/Debug-x64"
    RUNTIME_OUTPUT_DIRECTORY "/media/charlie/Data/dev/C++/SCESS/bin/linux/Debug-x64"
  )
endif()
target_include_directories("SCESS" PRIVATE
  $<$<CONFIG:Debug>:/media/charlie/Data/dev/C++/SCESS/src>
  $<$<CONFIG:Debug>:/usr/include>
  $<$<CONFIG:Debug>:/usr/local/include>
  $<$<CONFIG:Debug>:/media/charlie/Data/dev/Libraries/C++/json-develop/include>
  $<$<CONFIG:Debug>:/media/charlie/Data/dev/Libraries/C++/eigen-3.4.0>
  $<$<CONFIG:Debug>:/media/charlie/Data/dev/Libraries/C++/glm>
  $<$<CONFIG:Debug>:/media/charlie/Data/dev/Libraries/C++/ImGui>
)
target_compile_definitions("SCESS" PRIVATE
  $<$<CONFIG:Debug>:UNI_OPENGL>
  $<$<CONFIG:Debug>:_CONSOLE>
  $<$<CONFIG:Debug>:GLEW_STATIC>
  $<$<CONFIG:Debug>:CONFIG_DEBUG>
  $<$<CONFIG:Debug>:_DEBUG>
)
target_link_directories("SCESS" PRIVATE
)
find_package(PkgConfig REQUIRED)
pkg_check_modules(GTK REQUIRED IMPORTED_TARGET "gtk+-3.0")

target_link_libraries("SCESS"
  $<$<CONFIG:Debug>:libglfw3.a>
  $<$<CONFIG:Debug>:/usr/local/lib/libImGui.a>
  $<$<CONFIG:Debug>:/usr/local/lib/libboost_filesystem.a>
  $<$<CONFIG:Debug>:/usr/lib64/libGLEW.a>
  $<$<CONFIG:Debug>:GL>
  $<$<CONFIG:Debug>:m>
  $<$<CONFIG:Debug>:Xrandr>
  $<$<CONFIG:Debug>:Xi>
  $<$<CONFIG:Debug>:X11>
  $<$<CONFIG:Debug>:Xxf86vm>
  $<$<CONFIG:Debug>:pthread>
    PkgConfig::GTK
  $<$<CONFIG:Debug>:nfd.a>
  $<$<CONFIG:Debug>:/opt/etsf/lib/libxc.a>

)
#if(CMAKE_BUILD_TYPE STREQUAL Debug)
#  set_target_properties("SCESS" PROPERTIES COMPILE_FLAGS "$(pkg-config --cflags --libs gtk+-3.0) ")
#endif()
target_compile_options("SCESS" PRIVATE
  $<$<AND:$<CONFIG:Debug>,$<COMPILE_LANGUAGE:C>>:-m64>
  $<$<AND:$<CONFIG:Debug>,$<COMPILE_LANGUAGE:C>>:-g>
  $<$<AND:$<CONFIG:Debug>,$<COMPILE_LANGUAGE:CXX>>:-m64>
  $<$<AND:$<CONFIG:Debug>,$<COMPILE_LANGUAGE:CXX>>:-g>
  $<$<AND:$<CONFIG:Debug>,$<COMPILE_LANGUAGE:CXX>>:-std=c++20>
)
if(CMAKE_BUILD_TYPE STREQUAL Debug)
  set_target_properties("SCESS" PROPERTIES
    CXX_STANDARD 20
    CXX_STANDARD_REQUIRED YES
    CXX_EXTENSIONS NO
    POSITION_INDEPENDENT_CODE False
    INTERPROCEDURAL_OPTIMIZATION False
  )
endif()
if(CMAKE_BUILD_TYPE STREQUAL Release)
  set_target_properties("SCESS" PROPERTIES
    OUTPUT_NAME "SCESS"
    ARCHIVE_OUTPUT_DIRECTORY "/media/charlie/Data/dev/C++/SCESS/bin/linux/Release-x64"
    LIBRARY_OUTPUT_DIRECTORY "/media/charlie/Data/dev/C++/SCESS/bin/linux/Release-x64"
    RUNTIME_OUTPUT_DIRECTORY "/media/charlie/Data/dev/C++/SCESS/bin/linux/Release-x64"
  )
endif()
target_include_directories("SCESS" PRIVATE
  $<$<CONFIG:Release>:/media/charlie/Data/dev/C++/SCESS/src>
  $<$<CONFIG:Release>:/usr/include>
  $<$<CONFIG:Release>:/usr/local/include>
  $<$<CONFIG:Release>:/media/charlie/Data/dev/Libraries/C++/json-develop/include>
  $<$<CONFIG:Release>:/media/charlie/Data/dev/Libraries/C++/eigen-3.4.0>
  $<$<CONFIG:Release>:/media/charlie/Data/dev/Libraries/C++/glm>
  $<$<CONFIG:Release>:/media/charlie/Data/dev/Libraries/C++/ImGui>
)
target_compile_definitions("SCESS" PRIVATE
  $<$<CONFIG:Release>:UNI_OPENGL>
  $<$<CONFIG:Release>:_CONSOLE>
  $<$<CONFIG:Release>:GLEW_STATIC>
  $<$<CONFIG:Release>:CONFIG_RELEASE>
)
target_link_directories("SCESS" PRIVATE
)
target_link_libraries("SCESS"
  $<$<CONFIG:Release>:libglfw3.a>
  $<$<CONFIG:Release>:/usr/local/lib/libImGui.a>
  $<$<CONFIG:Release>:/usr/local/lib/libboost_filesystem.a>
  $<$<CONFIG:Release>:/usr/lib64/libGLEW.a>
  $<$<CONFIG:Release>:/usr/lib64/libGLEW.so>
  $<$<CONFIG:Release>:GL>
  $<$<CONFIG:Release>:m>
  $<$<CONFIG:Release>:Xrandr>
  $<$<CONFIG:Release>:Xi>
  $<$<CONFIG:Release>:X11>
  $<$<CONFIG:Release>:Xxf86vm>
  $<$<CONFIG:Release>:pthread>
    PkgConfig::GTK
  $<$<CONFIG:Release>:nfd.a>
  $<$<CONFIG:Release>:/opt/etsf/lib/libxc.a>

)
if(CMAKE_BUILD_TYPE STREQUAL Release)
  set_target_properties("SCESS" PROPERTIES COMPILE_FLAGS "")
endif()
target_compile_options("SCESS" PRIVATE
  $<$<AND:$<CONFIG:Release>,$<COMPILE_LANGUAGE:C>>:-m64>
  $<$<AND:$<CONFIG:Release>,$<COMPILE_LANGUAGE:C>>:-O2>
  $<$<AND:$<CONFIG:Release>,$<COMPILE_LANGUAGE:CXX>>:-m64>
  $<$<AND:$<CONFIG:Release>,$<COMPILE_LANGUAGE:CXX>>:-O2>
  $<$<AND:$<CONFIG:Release>,$<COMPILE_LANGUAGE:CXX>>:-std=c++20>
)
if(CMAKE_BUILD_TYPE STREQUAL Release)
  set_target_properties("SCESS" PROPERTIES
    CXX_STANDARD 20
    CXX_STANDARD_REQUIRED YES
    CXX_EXTENSIONS NO
    POSITION_INDEPENDENT_CODE False
    INTERPROCEDURAL_OPTIMIZATION False
  )
endif()
