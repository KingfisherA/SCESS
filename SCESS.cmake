add_executable("SCESS"
      "/opt/etsf/include"
      "src/SCESS/Application.cpp"
      "src/SCESS/Tests.cpp"
      "src/SCESS/Tests.h"
      "src/Utilities/Intrinsics.h"
      "src/Utilities/Resources.h"
      "src/Utilities/Utilities.cpp"
      "src/Utilities/Utilities.h"
)
if(CMAKE_BUILD_TYPE STREQUAL Debug)
  add_dependencies("SCESS"
    "QuantumSCESS"
  )
  set_target_properties("SCESS" PROPERTIES
    OUTPUT_NAME "SCESS"
    ARCHIVE_OUTPUT_DIRECTORY "/media/charlie/Data/dev/C++/SCESS/bin/linux/Debug-x64/SCESS"
    LIBRARY_OUTPUT_DIRECTORY "/media/charlie/Data/dev/C++/SCESS/bin/linux/Debug-x64/SCESS"
    RUNTIME_OUTPUT_DIRECTORY "/media/charlie/Data/dev/C++/SCESS/bin/linux/Debug-x64/SCESS"
  )
endif()
target_include_directories("SCESS" PRIVATE
  $<$<CONFIG:Debug>:/media/charlie/Data/dev/C++/SCESS/src>
  $<$<CONFIG:Debug>:/usr/include>
  $<$<CONFIG:Debug>:/usr/local/include>
  $<$<CONFIG:Debug>:/opt/etsf/include>
  $<$<CONFIG:Debug>:/media/charlie/Data/dev/Libraries/eigen-3.4.90>
  $<$<CONFIG:Debug>:/media/charlie/Data/dev/C++/GraphicsEngine/src>
)
target_compile_definitions("SCESS" PRIVATE
  $<$<CONFIG:Debug>:OS_LINUX>
  $<$<CONFIG:Debug>:CONFIG_DEBUG>
  $<$<CONFIG:Debug>:_DEBUG>
)
target_link_directories("SCESS" PRIVATE
  $<$<CONFIG:Debug>:/opt/etsf/lib>
  $<$<CONFIG:Debug>:/media/charlie/Data/dev/C++/GraphicsEngine/bin/linux/Release-x64/GraphicsEngine>
)
target_link_libraries("SCESS"
  $<$<CONFIG:Debug>:QuantumSCESS>
)
target_compile_options("SCESS" PRIVATE
  $<$<AND:$<CONFIG:Debug>,$<COMPILE_LANGUAGE:C>>:-m64>
  $<$<AND:$<CONFIG:Debug>,$<COMPILE_LANGUAGE:C>>:-O0>
  $<$<AND:$<CONFIG:Debug>,$<COMPILE_LANGUAGE:C>>:-g>
  $<$<AND:$<CONFIG:Debug>,$<COMPILE_LANGUAGE:CXX>>:-m64>
  $<$<AND:$<CONFIG:Debug>,$<COMPILE_LANGUAGE:CXX>>:-O0>
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
if(CMAKE_BUILD_TYPE STREQUAL Debug)
target_precompile_headers("SCESS" PUBLIC src/pch.h)
endif()
if(CMAKE_BUILD_TYPE STREQUAL Release)
  add_dependencies("SCESS"
    "QuantumSCESS"
  )
  set_target_properties("SCESS" PROPERTIES
    OUTPUT_NAME "SCESS"
    ARCHIVE_OUTPUT_DIRECTORY "/media/charlie/Data/dev/C++/SCESS/bin/linux/Release-x64/SCESS"
    LIBRARY_OUTPUT_DIRECTORY "/media/charlie/Data/dev/C++/SCESS/bin/linux/Release-x64/SCESS"
    RUNTIME_OUTPUT_DIRECTORY "/media/charlie/Data/dev/C++/SCESS/bin/linux/Release-x64/SCESS"
  )
endif()
target_include_directories("SCESS" PRIVATE
  $<$<CONFIG:Release>:/media/charlie/Data/dev/C++/SCESS/src>
  $<$<CONFIG:Release>:/usr/include>
  $<$<CONFIG:Release>:/usr/local/include>
  $<$<CONFIG:Release>:/opt/etsf/include>
  $<$<CONFIG:Release>:/media/charlie/Data/dev/Libraries/eigen-3.4.90>
  $<$<CONFIG:Release>:/media/charlie/Data/dev/C++/GraphicsEngine/src>
)
target_compile_definitions("SCESS" PRIVATE
  $<$<CONFIG:Release>:OS_LINUX>
  $<$<CONFIG:Release>:CONFIG_RELEASE>
)
target_link_directories("SCESS" PRIVATE
  $<$<CONFIG:Release>:/opt/etsf/lib>
  $<$<CONFIG:Release>:/media/charlie/Data/dev/C++/GraphicsEngine/bin/linux/Release-x64/GraphicsEngine>
)
target_link_libraries("SCESS"
  $<$<CONFIG:Release>:QuantumSCESS>
)
target_compile_options("SCESS" PRIVATE
  $<$<AND:$<CONFIG:Release>,$<COMPILE_LANGUAGE:C>>:-m64>
  $<$<AND:$<CONFIG:Release>,$<COMPILE_LANGUAGE:C>>:-O3>
  $<$<AND:$<CONFIG:Release>,$<COMPILE_LANGUAGE:C>>:-mavx2>
  $<$<AND:$<CONFIG:Release>,$<COMPILE_LANGUAGE:CXX>>:-m64>
  $<$<AND:$<CONFIG:Release>,$<COMPILE_LANGUAGE:CXX>>:-O3>
  $<$<AND:$<CONFIG:Release>,$<COMPILE_LANGUAGE:CXX>>:-mavx2>
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
if(CMAKE_BUILD_TYPE STREQUAL Release)
target_precompile_headers("SCESS" PUBLIC src/pch.h)
endif()