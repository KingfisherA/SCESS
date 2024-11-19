cd ..
(cd ./Libraries && ./Build_Libraries.sh)

#cmake -DCMAKE_BUILD_TYPE=Debug -DCMAKE_MAKE_PROGRAM=ninja -DCMAKE_C_COMPILER=clang -DCMAKE_CXX_COMPILER=clang++ -G Ninja -S ./ -B ./cmake-build-debug
#cmake --build ./cmake-build-debug --target all -j 10

cmake -DCMAKE_BUILD_TYPE=Release -DCMAKE_MAKE_PROGRAM=ninja -DCMAKE_C_COMPILER=clang -DCMAKE_CXX_COMPILER=clang++ -G Ninja -S ./ -B ./cmake-build-release
cmake --build ./cmake-build-release --target all -j 10

