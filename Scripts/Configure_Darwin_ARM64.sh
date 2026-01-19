#!/bin/bash
while getopts ":c:" opt; do
  case $opt in
    c) CONFIG="$OPTARG"
    ;;
    \?) echo "Invalid option -$OPTARG" >&2
    exit 1
    ;;
  esac

  case $OPTARG in
    -*) echo "Option $opt needs a valid argument"
    exit 1
    ;;
  esac
done

cmake -DCMAKE_BUILD_TYPE=$CONFIG -DCMAKE_MAKE_PROGRAM=ninja -DCMAKE_C_COMPILER=clang -DCMAKE_CXX_COMPILER=clang++ -G Ninja -S ../ -B ../cmake-build-"${CONFIG}"-darwin

