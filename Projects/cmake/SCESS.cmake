add_executable("SCESS"
  "${SCESS_SRC_DIR}/SCESS/Application.cpp"
  "${SCESS_SRC_DIR}/SCESS/Tests.cpp"
  "${SCESS_SRC_DIR}/SCESS/Tests.h"
)

add_dependencies("SCESS"
    "QuantumSCESS"
)
target_compile_options("SCESS" PRIVATE
  -m64
  -std=c++20
)
target_include_directories("SCESS" PRIVATE
  ${SCESS_SRC_DIR}
  ${SCESS_LIB_DIR}/GraphicsEngine-1.0.0/include
  ${SCESS_LIB_DIR}/cppitertools
  ${SCESS_LIB_DIR}/JSON/include
  ${SCESS_LIB_DIR}/Eigen
)
target_compile_definitions("SCESS" PRIVATE
  OS_LINUX
  $<$<CONFIG:DEBUG>:CONFIG_DEBUG>
  $<$<CONFIG:DEBUG>:_DEBUG>
  $<$<CONFIG:RELEASE>:CONFIG_RELEASE>
)
target_link_libraries("SCESS" PRIVATE
  QuantumSCESS
  ${LIB_GraphicsEngine}
  Libint2::libint2
)
target_precompile_headers("SCESS" PUBLIC ./src/pch.h)

if(CMAKE_BUILD_TYPE STREQUAL Debug)
  set_target_properties("SCESS" PROPERTIES
    CXX_STANDARD 20
    CXX_STANDARD_REQUIRED YES
    CXX_EXTENSIONS NO
    POSITION_INDEPENDENT_CODE False
    INTERPROCEDURAL_OPTIMIZATION False
    OUTPUT_NAME "SCESS"
    ARCHIVE_OUTPUT_DIRECTORY "${SCESS_BIN_DIR}/Debug-x64"
    LIBRARY_OUTPUT_DIRECTORY "${SCESS_BIN_DIR}/Debug-x64"
    RUNTIME_OUTPUT_DIRECTORY "${SCESS_BIN_DIR}/Debug-x64"
  )
  target_compile_options("SCESS" PRIVATE
    -O0
    -g
  ) 
endif()

if(CMAKE_BUILD_TYPE STREQUAL Release)
  set_target_properties("SCESS" PROPERTIES
    CXX_STANDARD 20
    CXX_STANDARD_REQUIRED YES
    CXX_EXTENSIONS NO
    POSITION_INDEPENDENT_CODE False
    INTERPROCEDURAL_OPTIMIZATION False
    OUTPUT_NAME "SCESS"
    ARCHIVE_OUTPUT_DIRECTORY "${SCESS_BIN_DIR}/Release-x64"
    LIBRARY_OUTPUT_DIRECTORY "${SCESS_BIN_DIR}/Release-x64"
    RUNTIME_OUTPUT_DIRECTORY "${SCESS_BIN_DIR}/Release-x64"
  )
  target_compile_options("SCESS" PRIVATE
    -O2
    -mavx2
  ) 
endif()