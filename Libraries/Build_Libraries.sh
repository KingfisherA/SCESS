if [ ! -d "./LibXC/build" ]; then
	cmake -H./LibXC/ -B./LibXC/build
fi
(cd ./LibXC/build && make)
if [ ! -d "./LibXC/lib" ]; then
	mkdir "./LibXC/lib"
	mkdir "./LibXC/lib/include"
fi
find "./LibXC/build/" -name \*.a -exec cp {} "./LibXC/lib" \;
find "./LibXC/build/" -name \*.h -exec cp {} "./LibXC/lib/include" \;
find "./LibXC/build/src" -name \*.h -exec cp {} "./LibXC/lib/include" \;

(cd ./nfd/build/gmake_linux && make config=release_x64)

mkdir -p -- ../bin/linux/Debug-x64/
mkdir -p  -- ../bin/linux/Release-x64/

cp --update=none ./imgui.ini ../bin/linux/Release-x64/imgui.ini
cp --update=none ./imgui.ini ../bin/linux/Debug-x64/imgui.ini
cp --update=none ./settings.json ../bin/linux/Release-x64/settings.json
cp --update=none ./settings.json ../bin/linux/Debug-x64/settings.json

cp --update=none ./GraphicsEngine-1.0.0/libGraphicsEngine.so ../bin/linux/Release-x64/libGraphicsEngine.so
cp --update=none ./GraphicsEngine-1.0.0/libGraphicsEngine.so ../bin/linux/Debug-x64/libGraphicsEngine.so
