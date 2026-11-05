# Settings shared by the SCESS targets.
#
#   scess_target_defaults(<target>)
#
# C++20 without extensions, the shared MSVC runtime, the CONFIG_* / OS_* defines used by the code,
# and an output directory of build/<target> (e.g. build/ViewerSCESS/ViewerSCESS). Optimisation and
# debug-info flags come from CMake's per-configuration defaults (Debug: -g, Release: -O3 -DNDEBUG).
function(scess_target_defaults TARGET)
  set_target_properties(${TARGET} PROPERTIES
    CXX_STANDARD 20
    CXX_STANDARD_REQUIRED ON
    CXX_EXTENSIONS OFF
    MSVC_RUNTIME_LIBRARY "MultiThreaded$<$<CONFIG:Debug>:Debug>DLL"
    ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/${TARGET}
    LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/${TARGET}
    RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/${TARGET}
  )
  target_compile_definitions(${TARGET} PRIVATE
    LIBINT2_CPLUSPLUS_STD=2020
    $<$<PLATFORM_ID:Linux>:OS_LINUX>
    $<$<PLATFORM_ID:Windows>:OS_WINDOWS>
    $<$<CONFIG:Debug>:CONFIG_DEBUG _DEBUG>
    $<$<CONFIG:Release>:CONFIG_RELEASE>
  )
endfunction()
