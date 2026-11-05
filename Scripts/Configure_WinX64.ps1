# Configures SCESS into .\build on Windows (dependencies from C:\Packages and vcpkg)
#   Scripts\Configure_WinX64.ps1 -config Release
param ($config = "Release")
$root = Split-Path -Parent $PSScriptRoot
cmake -DCMAKE_PREFIX_PATH="C:/Packages" -DCMAKE_BUILD_TYPE="$config" -DCMAKE_TOOLCHAIN_FILE="${env:VCPKG_ROOT}/scripts/buildsystems/vcpkg.cmake" -DVCPKG_TARGET_TRIPLET=x64-windows -S "$root" -B "$root/build"
