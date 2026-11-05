#!/bin/sh
# Configures SCESS into ./build and builds it (from any directory)
#
#   Scripts/Build.sh [Debug|Release] [cmake --build options...]
#
# e.g. Scripts/Build.sh Release --target ViewerSCESSApp
set -eu

ROOT="$(cd "$(dirname "$0")/.." && pwd)"
CONFIG="${1:-Release}"
[ $# -gt 0 ] && shift

# cmake.sh adds Homebrew's LAPACK to CMAKE_PREFIX_PATH and uses ccache when available
"$ROOT/cmake.sh" -S "$ROOT" -B "$ROOT/build" -G Ninja -DCMAKE_BUILD_TYPE="$CONFIG"
"$ROOT/cmake.sh" --build "$ROOT/build" "$@"
