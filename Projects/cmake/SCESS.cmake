add_executable("SCESS"
  "src/SCESS/Application.cpp"
  "src/SCESS/Tests.cpp"
  "src/SCESS/Tests.h"
)
if(CMAKE_BUILD_TYPE STREQUAL Debug)
  add_dependencies("SCESS"
    "QuantumSCESS"
  )
  set_target_properties("SCESS" PROPERTIES
    OUTPUT_NAME "SCESS"
    ARCHIVE_OUTPUT_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/bin/linux/Debug-x64"
    LIBRARY_OUTPUT_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/bin/linux/Debug-x64"
    RUNTIME_OUTPUT_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/bin/linux/Debug-x64"
  )
endif()
target_include_directories("SCESS" PRIVATE
  $<$<CONFIG:Debug>:${CMAKE_CURRENT_SOURCE_DIR}/src>
  $<$<CONFIG:Debug>:${CMAKE_CURRENT_SOURCE_DIR}/Libraries/GraphicsEngine-1.0.0/include>
  $<$<CONFIG:Debug>:${CMAKE_CURRENT_SOURCE_DIR}/Libraries/glm>
  $<$<CONFIG:Debug>:${CMAKE_CURRENT_SOURCE_DIR}/Libraries/cppitertools>
  $<$<CONFIG:Debug>:${CMAKE_CURRENT_SOURCE_DIR}/Libraries/JSON/include>
  $<$<CONFIG:Debug>:${CMAKE_CURRENT_SOURCE_DIR}/Libraries/Eigen>
  $<$<CONFIG:Debug>:${CMAKE_CURRENT_SOURCE_DIR}/Libraries/LibInt2/include>
  $<$<CONFIG:Debug>:${CMAKE_CURRENT_SOURCE_DIR}/Libraries/LibXC/lib/include>
)
target_compile_definitions("SCESS" PRIVATE
  $<$<CONFIG:Debug>:OS_LINUX>
  $<$<CONFIG:Debug>:CONFIG_DEBUG>
  $<$<CONFIG:Debug>:_DEBUG>
)
target_link_directories("SCESS" PRIVATE
  $<$<CONFIG:Debug>:${CMAKE_CURRENT_SOURCE_DIR}/Libraries/LibXC/lib>
  $<$<CONFIG:Debug>:${CMAKE_CURRENT_SOURCE_DIR}/Libraries/GraphicsEngine-1.0.0>
)
target_link_libraries("SCESS"
  $<$<CONFIG:Debug>:QuantumSCESS>
  $<$<CONFIG:Debug>:xc.a>
  $<$<CONFIG:Debug>:GraphicsEngine>
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
target_precompile_headers("SCESS" PUBLIC ./src/pch.h)
endif()
if(CMAKE_BUILD_TYPE STREQUAL Release)
  add_dependencies("SCESS"
    "QuantumSCESS"
  )
  set_target_properties("SCESS" PROPERTIES
    OUTPUT_NAME "SCESS"
    ARCHIVE_OUTPUT_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/bin/linux/Release-x64"
    LIBRARY_OUTPUT_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/bin/linux/Release-x64"
    RUNTIME_OUTPUT_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/bin/linux/Release-x64"
  )
endif()
target_include_directories("SCESS" PRIVATE
  $<$<CONFIG:Release>:${CMAKE_CURRENT_SOURCE_DIR}/src>
  $<$<CONFIG:Release>:${CMAKE_CURRENT_SOURCE_DIR}/Libraries/GraphicsEngine-1.0.0/include>
  $<$<CONFIG:Release>:${CMAKE_CURRENT_SOURCE_DIR}/Libraries/glm>
  $<$<CONFIG:Release>:${CMAKE_CURRENT_SOURCE_DIR}/Libraries/cppitertools>
  $<$<CONFIG:Release>:${CMAKE_CURRENT_SOURCE_DIR}/Libraries/JSON/include>
  $<$<CONFIG:Release>:${CMAKE_CURRENT_SOURCE_DIR}/Libraries/Eigen>
  $<$<CONFIG:Release>:${CMAKE_CURRENT_SOURCE_DIR}/Libraries/LibInt2/include>
  $<$<CONFIG:Release>:${CMAKE_CURRENT_SOURCE_DIR}/Libraries/LibXC/lib/include>
)
target_compile_definitions("SCESS" PRIVATE
  $<$<CONFIG:Release>:OS_LINUX>
  $<$<CONFIG:Release>:CONFIG_RELEASE>
)
target_link_directories("SCESS" PRIVATE
  $<$<CONFIG:Release>:${CMAKE_CURRENT_SOURCE_DIR}/Libraries/LibXC/lib>
  $<$<CONFIG:Release>:${CMAKE_CURRENT_SOURCE_DIR}/Libraries/GraphicsEngine-1.0.0>
)
target_link_libraries("SCESS"
  $<$<CONFIG:Release>:QuantumSCESS>
  $<$<CONFIG:Release>:xc.a>
  $<$<CONFIG:Release>:GraphicsEngine>
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
target_precompile_headers("SCESS" PUBLIC ./src/pch.h)
endif()