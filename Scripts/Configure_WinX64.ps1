param ($config)
cmake -DCMAKE_PREFIX_PATH="C:/Packages" -DEIGEN3_INCLUDE_DIR="./Libraries/Eigen" -DCMAKE_BUILD_TYPE="$($config)" -DCMAKE_TOOLCHAIN_FILE="${env:VCPKG_ROOT}/scripts/buildsystems/vcpkg.cmake" -DVCPKG_TARGET_TRIPLET=x64-windows -S ../ -B ../cmake-build-$($config.ToLower())-win
