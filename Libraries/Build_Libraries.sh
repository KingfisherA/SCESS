GNU_COPY="/opt/homebrew/Cellar/coreutils/9.7/bin/gcp"

if [ ! -d "./LibXC/build" ]; then
	cmake -H./LibXC/ -B./LibXC/build
fi
(cd ./LibXC/build && make -j 10)
if [ ! -d "./LibXC/lib" ]; then
	mkdir "./LibXC/lib"
	mkdir "./LibXC/lib/include"
fi
find "./LibXC/build/" -name \*.a -exec $GNU_COPY {} "./LibXC/lib" \;
find "./LibXC/build/" -name \*.h -exec $GNU_COPY {} "./LibXC/lib/include" \;
find "./LibXC/build/src" -name \*.h -exec $GNU_COPY {} "./LibXC/lib/include" \;

(cd ./nfd/build/gmake_linux && make config=release_x64)

cmake -DIMGUI_INCLUDE_DIRS="/Users/charlie/dev/C++/GraphicsEngine/include/ImGui/" -DCMAKE_BUILD_TYPE=Release -DCMAKE_MAKE_PROGRAM=ninja -DCMAKE_C_COMPILER=clang -DCMAKE_CXX_COMPILER=clang++ -G Ninja -S./ImPlot_Build/ -B./ImPlot/cmake-build-release
cmake --build ./ImPlot/cmake-build-release --target all -j 10  


# mkdir -p -- ../bin/linux/Debug-x64/
# mkdir -p  -- ../bin/linux/Release-x64/

# $GNU_COPY --update=none ./imgui.ini ../bin/linux/Release-x64/imgui.ini
# $GNU_COPY --update=none ./imgui.ini ../bin/linux/Debug-x64/imgui.ini
# $GNU_COPY --update=none ./settings.json ../bin/linux/Release-x64/settings.json
# $GNU_COPY --update=none ./settings.json ../bin/linux/Debug-x64/settings.json

# $GNU_COPY --update=none ./GraphicsEngine-1.0.0/libGraphicsEngine.so ../bin/linux/Release-x64/libGraphicsEngine.so
# $GNU_COPY --update=none ./GraphicsEngine-1.0.0/libGraphicsEngine.so ../bin/linux/Debug-x64/libGraphicsEngine.so
