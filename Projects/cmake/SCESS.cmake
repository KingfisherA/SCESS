add_executable("SCESS"
  "${SCESS_SRC_DIR}/Resources/AssetFolder.cpp"
  "${SCESS_SRC_DIR}/Resources/AssetFolder.h"
  "${SCESS_SRC_DIR}/SCESS/Application.cpp"
  "${SCESS_SRC_DIR}/SCESS/Tests.cpp"
  "${SCESS_SRC_DIR}/SCESS/Tests.h"
        "${SCESS_SRC_DIR}/SCESS/TestEinsum.cpp"
)

add_dependencies("SCESS"
    "QuantumSCESS"
)
target_include_directories("SCESS" PRIVATE
  ${SCESS_SRC_DIR}
  ${GENG_INCLUDE_DIRS}
  ${IMGUI_INCLUDE_DIRS}
  ${SCESS_LIB_DIR}/cppitertools
  ${SCESS_LIB_DIR}/JSON/include
  ${SCESS_LIB_DIR}/Eigen
)
target_compile_definitions("SCESS" PRIVATE
  LIBINT2_CPLUSPLUS_STD=2020
  $<$<PLATFORM_ID:Linux>:OS_LINUX>
  $<$<PLATFORM_ID:Windows>:OS_WINDOWS>
  $<$<CONFIG:DEBUG>:CONFIG_DEBUG>
  $<$<CONFIG:DEBUG>:_DEBUG>
  $<$<CONFIG:RELEASE>:CONFIG_RELEASE>
)
target_link_libraries("SCESS" PRIVATE
  fmt
  QuantumSCESS
  ${GENG_TOOLS_LIBRARIES}
  Boost::filesystem
  $<TARGET_NAME_IF_EXISTS:Libint2::int2>
  $<TARGET_NAME_IF_EXISTS:Libint2::libint2>
)

#target_precompile_headers("SCESS" PUBLIC ./src/pch.h)

set(COMPILE_OPTIONS)
if((CMAKE_CXX_COMPILER_ID STREQUAL Clang) OR 
   (CMAKE_CXX_COMPILER_ID STREQUAL GNU))
    list(APPEND COMPILE_OPTIONS
      -m64
      -fPIC
      -std=c++20
    )
    if(CMAKE_BUILD_TYPE STREQUAL Debug)
      list(APPEND COMPILE_OPTIONS
        -g
        -O0
      )
    elseif(CMAKE_BUILD_TYPE STREQUAL Release)
      list(APPEND COMPILE_OPTIONS
        -O3
        -mavx2
      )
    endif()
elseif(CMAKE_CXX_COMPILER_ID STREQUAL MSVC)
  list(APPEND COMPILE_OPTIONS
      /std:c++20
    )
    if(CMAKE_BUILD_TYPE STREQUAL Debug)
      list(APPEND COMPILE_OPTIONS
        /Zi
        /Od
      )
    elseif(CMAKE_BUILD_TYPE STREQUAL Release)
      list(APPEND COMPILE_OPTIONS /O2)
    endif()
endif()
target_compile_options("SCESS" PRIVATE ${COMPILE_OPTIONS})

set_target_properties("SCESS" PROPERTIES
    CXX_STANDARD 20
    CXX_STANDARD_REQUIRED YES
    CXX_EXTENSIONS NO
    POSITION_INDEPENDENT_CODE True
    INTERPROCEDURAL_OPTIMIZATION False
    MSVC_RUNTIME_LIBRARY MultiThreaded$<$<CONFIG:Debug>:Debug>DLL
    OUTPUT_NAME "SCESS")

if(CMAKE_CONFIGURATION_TYPES)
  set(DEBUG_BIN_DIR "${SCESS_BIN_DIR}")
  set(RELEASE_BIN_DIR "${SCESS_BIN_DIR}")
else()
  set(DEBUG_BIN_DIR "${SCESS_BIN_DIR}/Debug")
  set(RELEASE_BIN_DIR "${SCESS_BIN_DIR}/Release")
endif()

if(CMAKE_BUILD_TYPE STREQUAL Debug)
  set_target_properties("SCESS" PROPERTIES
    ARCHIVE_OUTPUT_DIRECTORY "${DEBUG_BIN_DIR}"
    LIBRARY_OUTPUT_DIRECTORY "${DEBUG_BIN_DIR}"
    RUNTIME_OUTPUT_DIRECTORY "${DEBUG_BIN_DIR}"
  )
elseif(CMAKE_BUILD_TYPE STREQUAL Release)
  set_target_properties("SCESS" PROPERTIES
    ARCHIVE_OUTPUT_DIRECTORY "${RELEASE_BIN_DIR}"
    LIBRARY_OUTPUT_DIRECTORY "${RELEASE_BIN_DIR}"
    RUNTIME_OUTPUT_DIRECTORY "${RELEASE_BIN_DIR}"
  )
endif()
