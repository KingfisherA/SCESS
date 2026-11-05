#!/bin/sh
# Homebrew's LAPACK (keg-only), then Homebrew itself so its packages win over others found via PATH
# (e.g. a conda fmt whose version doesn't match the Homebrew headers on the include path)
export CMAKE_PREFIX_PATH=/opt/homebrew/opt/lapack:/opt/homebrew

# Route compiles through ccache when it is installed; CMake reads these at configure time
if command -v ccache >/dev/null 2>&1; then
    export CMAKE_C_COMPILER_LAUNCHER=ccache
    export CMAKE_CXX_COMPILER_LAUNCHER=ccache
fi

exec cmake "$@"
