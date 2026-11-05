#!/bin/sh
# Installs the SCESS (and Dahlia) build dependencies on Debian/Ubuntu. Libint2 is not packaged with
# C++ support and needs building from source: https://github.com/evaleev/libint
apt install cmake ninja-build clang glslc \
    libboost-all-dev libblas-dev liblapack-dev liblapacke-dev libfmt-dev libxc-dev \
    libglfw3-dev libglew-dev libglm-dev libfreetype-dev libgtk-3-dev
