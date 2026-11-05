# ViewerSCESS: the graphical viewer, built into build/ViewerSCESS
#
#   ViewerSCESS      plain executable (default build)
#   ViewerSCESSApp   macOS only: ViewerSCESS.app with the SCESS icon, the Assets folder and every
#                    non-system library embedded (cmake --build build --target ViewerSCESSApp)

# The viewer's sources are compiled once into an object library that both executables link
add_library(ViewerSCESSObjects OBJECT
  ${SCESS_SRC_DIR}/Resources/ShaderSources.cpp
  ${SCESS_SRC_DIR}/Resources/ShaderSources.h
  ${SCESS_SRC_DIR}/Utilities/Assets.h
  ${SCESS_SRC_DIR}/Utilities/Intrinsics.h
  ${SCESS_SRC_DIR}/Utilities/Utilities.cpp
  ${SCESS_SRC_DIR}/Utilities/Utilities.h
  ${SCESS_SRC_DIR}/Viewer/Application/Application.cpp
  ${SCESS_SRC_DIR}/Viewer/Application/Application.h
  ${SCESS_SRC_DIR}/Viewer/GUI/AtomDraw.cpp
  ${SCESS_SRC_DIR}/Viewer/GUI/AtomDraw.h
  ${SCESS_SRC_DIR}/Viewer/GUI/CalculationLayer.cpp
  ${SCESS_SRC_DIR}/Viewer/GUI/CalculationLayer.h
  ${SCESS_SRC_DIR}/Viewer/GUI/ImSpinner.h
  ${SCESS_SRC_DIR}/Viewer/Isosurfaces/IsosurfaceConstants.cpp
  ${SCESS_SRC_DIR}/Viewer/Isosurfaces/Isosurfaces.cpp
  ${SCESS_SRC_DIR}/Viewer/Isosurfaces/Isosurfaces.h
)
scess_target_defaults(ViewerSCESSObjects)
target_link_libraries(ViewerSCESSObjects PUBLIC
  QuantumSCESS
  Dahlia
  DahliaTools
  ImPlot
  nfd
)

# Vulkan-style GLSL compiled to SPIR-V and embedded as .inc files (needs glslc)
dahlia_compile_shaders(ViewerSCESSObjects
  ${SCESS_SRC_DIR}/Resources/Shaders/MoleculeDiffuse.vert
  ${SCESS_SRC_DIR}/Resources/Shaders/MoleculeDiffuse.frag
  ${SCESS_SRC_DIR}/Resources/Shaders/IsosurfaceDiffuse.vert
  ${SCESS_SRC_DIR}/Resources/Shaders/IsosurfaceDiffuse.frag
  ${SCESS_SRC_DIR}/Resources/Shaders/SolidColour.vert
  ${SCESS_SRC_DIR}/Resources/Shaders/SolidColour.frag
)

# ---- ViewerSCESS ----------------------------------------------------------------------------------
add_executable(ViewerSCESS)
scess_target_defaults(ViewerSCESS)
target_link_libraries(ViewerSCESS PRIVATE ViewerSCESSObjects)

if(WIN32)
  # Windows has no rpath: put the runtime DLLs next to the executable
  add_custom_command(TARGET ViewerSCESS POST_BUILD
    COMMAND ${CMAKE_COMMAND} -E copy_if_different
            $<TARGET_RUNTIME_DLLS:ViewerSCESS> $<TARGET_FILE_DIR:ViewerSCESS>
    COMMAND_EXPAND_LISTS
  )
endif()

# ---- ViewerSCESS.app ------------------------------------------------------------------------------
if(APPLE)
  add_executable(ViewerSCESSApp MACOSX_BUNDLE EXCLUDE_FROM_ALL)
  scess_target_defaults(ViewerSCESSApp)
  target_link_libraries(ViewerSCESSApp PRIVATE ViewerSCESSObjects)
  set_target_properties(ViewerSCESSApp PROPERTIES
    OUTPUT_NAME ViewerSCESS
    # Next to the plain executable: build/ViewerSCESS/ViewerSCESS.app
    RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/ViewerSCESS
    MACOSX_BUNDLE_INFO_PLIST ${PROJECT_SOURCE_DIR}/cmake/ViewerSCESS.plist.in
    MACOSX_BUNDLE_BUNDLE_NAME SCESS
    MACOSX_BUNDLE_GUI_IDENTIFIER com.kingfishera.SCESS
    MACOSX_BUNDLE_ICON_FILE "SCESS Icon.icns"
    MACOSX_BUNDLE_SHORT_VERSION_STRING ${PROJECT_VERSION}
    MACOSX_BUNDLE_BUNDLE_VERSION ${PROJECT_VERSION}
    # Load libraries from the bundle's Frameworks folder rather than the build tree
    BUILD_WITH_INSTALL_RPATH ON
    INSTALL_RPATH @executable_path/../Frameworks
    # Relink (re-running the POST_BUILD steps below) when the icon is regenerated
    LINK_DEPENDS "${SCESS_ASSETS_DIR}/Icon/SCESS Icon.icns"
  )

  # BundleDependencies.cmake resolves the libraries to embed from the plain executable
  add_dependencies(ViewerSCESSApp ViewerSCESS)

  set(SCESS_BUNDLE_RESOURCES $<TARGET_BUNDLE_CONTENT_DIR:ViewerSCESSApp>/Resources)
  add_custom_command(TARGET ViewerSCESSApp POST_BUILD
    # The icon is copied here rather than through MACOSX_PACKAGE_LOCATION, which Ninja only runs
    # before compiling the target's own sources (it has none). Regenerate it from the SVG with
    # Scripts/Generate_Icons.sh
    COMMAND ${CMAKE_COMMAND} -E copy_if_different "${SCESS_ASSETS_DIR}/Icon/SCESS Icon.icns"
            "${SCESS_BUNDLE_RESOURCES}/SCESS Icon.icns"
    # Assets are found at Contents/Resources/Assets (see Resources/AssetFolder.cpp)
    COMMAND ${CMAKE_COMMAND} -E rm -rf ${SCESS_BUNDLE_RESOURCES}/Assets
    COMMAND ${CMAKE_COMMAND} -E copy_directory ${SCESS_ASSETS_DIR} ${SCESS_BUNDLE_RESOURCES}/Assets
    # Embeds libQuantumSCESS, libDahlia and their non-system dependencies, then signs the bundle
    COMMAND ${CMAKE_COMMAND}
            -DAPP=$<TARGET_BUNDLE_DIR:ViewerSCESSApp>
            -DREFERENCE=$<TARGET_FILE:ViewerSCESS>
            -P ${PROJECT_SOURCE_DIR}/cmake/BundleDependencies.cmake
    VERBATIM
  )
endif()
