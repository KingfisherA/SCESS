# Makes a macOS .app self-contained: copies every non-system dylib the executable needs into
# Contents/Frameworks, points all references at @rpath/<name> (the executable's rpath is
# @executable_path/../Frameworks), then ad-hoc signs the bundle.
#
#   cmake -DAPP=<path/to/X.app> -DREFERENCE=<executable> -P BundleDependencies.cmake
#
# REFERENCE is the same program built outside the bundle: it still has the build tree's rpaths, so
# its @rpath dependencies (libQuantumSCESS, libDahlia, libint2) can be resolved.

if(NOT APP OR NOT IS_DIRECTORY "${APP}")
  message(FATAL_ERROR "BundleDependencies: APP must be an existing .app bundle (got '${APP}')")
endif()
if(NOT REFERENCE OR NOT EXISTS "${REFERENCE}")
  message(FATAL_ERROR "BundleDependencies: REFERENCE must be an existing executable (got '${REFERENCE}')")
endif()

get_filename_component(APP_NAME "${APP}" NAME_WE)
set(EXECUTABLE "${APP}/Contents/MacOS/${APP_NAME}")
set(FRAMEWORKS "${APP}/Contents/Frameworks")

file(REMOVE_RECURSE "${FRAMEWORKS}")
file(MAKE_DIRECTORY "${FRAMEWORKS}")

file(GET_RUNTIME_DEPENDENCIES
  EXECUTABLES "${REFERENCE}"
  RESOLVED_DEPENDENCIES_VAR DEPENDENCIES
  UNRESOLVED_DEPENDENCIES_VAR UNRESOLVED
  POST_EXCLUDE_REGEXES "^/usr/lib/" "^/System/"
)
if(UNRESOLVED)
  message(FATAL_ERROR "BundleDependencies: could not find ${UNRESOLVED}")
endif()

# Copy each library (following symlinks) under the name it is referenced by
set(BUNDLED_NAMES)
set(BUNDLED_FILES)
foreach(DEPENDENCY ${DEPENDENCIES})
  get_filename_component(NAME "${DEPENDENCY}" NAME)
  file(REAL_PATH "${DEPENDENCY}" SOURCE)
  file(COPY_FILE "${SOURCE}" "${FRAMEWORKS}/${NAME}")
  file(CHMOD "${FRAMEWORKS}/${NAME}" PERMISSIONS OWNER_READ OWNER_WRITE GROUP_READ WORLD_READ)
  list(APPEND BUNDLED_NAMES "${NAME}")
  list(APPEND BUNDLED_FILES "${FRAMEWORKS}/${NAME}")
endforeach()

function(run)
  execute_process(COMMAND ${ARGN} RESULT_VARIABLE RESULT OUTPUT_VARIABLE OUTPUT ERROR_VARIABLE OUTPUT)
  if(NOT RESULT EQUAL 0)
    string(REPLACE ";" " " COMMAND "${ARGN}")
    message(FATAL_ERROR "BundleDependencies: '${COMMAND}' failed:\n${OUTPUT}")
  endif()
endfunction()

# Rewrites one Mach-O file: bundled libraries are referenced as @rpath/<name>, and a library's own
# rpaths (absolute build/Homebrew directories) are replaced by @loader_path (Contents/Frameworks).
function(relink FILE IS_LIBRARY)
  execute_process(COMMAND otool -L "${FILE}" OUTPUT_VARIABLE LINKS)
  string(REPLACE "\n" ";" LINKS "${LINKS}")
  list(REMOVE_AT LINKS 0) # "<file>:"
  set(ARGS)
  foreach(LINE ${LINKS})
    if(LINE MATCHES "^[ \t]+([^ \t].*) \\(compatibility version")
      set(INSTALL_NAME "${CMAKE_MATCH_1}")
      get_filename_component(NAME "${INSTALL_NAME}" NAME)
      if(NAME IN_LIST BUNDLED_NAMES AND NOT INSTALL_NAME STREQUAL "@rpath/${NAME}")
        list(APPEND ARGS -change "${INSTALL_NAME}" "@rpath/${NAME}")
      endif()
    endif()
  endforeach()

  if(IS_LIBRARY)
    get_filename_component(NAME "${FILE}" NAME)
    list(APPEND ARGS -id "@rpath/${NAME}")
    execute_process(COMMAND otool -l "${FILE}" OUTPUT_VARIABLE LOAD_COMMANDS)
    string(REGEX MATCHALL "cmd LC_RPATH[^\n]*\n[^\n]*\n[ \t]+path [^\n]+ \\(offset" RPATHS "${LOAD_COMMANDS}")
    set(HAS_LOADER_PATH FALSE)
    foreach(RPATH ${RPATHS})
      string(REGEX REPLACE ".*path (.+) \\(offset" "\\1" RPATH "${RPATH}")
      if(RPATH STREQUAL "@loader_path")
        set(HAS_LOADER_PATH TRUE)
      else()
        list(APPEND ARGS -delete_rpath "${RPATH}")
      endif()
    endforeach()
    if(NOT HAS_LOADER_PATH)
      list(APPEND ARGS -add_rpath @loader_path)
    endif()
  endif()

  if(ARGS)
    # install_name_tool warns that it invalidates the signature; the bundle is re-signed below
    run(install_name_tool ${ARGS} "${FILE}")
  endif()
endfunction()

relink("${EXECUTABLE}" FALSE)
foreach(FILE ${BUNDLED_FILES})
  relink("${FILE}" TRUE)
endforeach()

# Ad-hoc sign the libraries, then the bundle, so its contents are sealed (and allowed to run on
# Apple Silicon after being modified above)
foreach(FILE ${BUNDLED_FILES})
  run(codesign --force --sign - "${FILE}")
endforeach()
run(codesign --force --sign - "${APP}")

list(LENGTH BUNDLED_FILES COUNT)
message(STATUS "Bundled ${COUNT} libraries into ${FRAMEWORKS}")
