cmake -B build -DCMAKE_TOOLCHAIN_FILE=C:/DevTools/vcpkg/scripts/buildsystem/vcpkg.cmake -DVCPKG_TARGET_TRIPLET=x64-windows
cmake --build build
./build/Debug/replace.exe --lef ./test/ispd18_test1.input.lef --def ./test/ispd18_test1.input.def

