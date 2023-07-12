cmake -DUSE_CIMG=ON -B build
cmake --build build
./build/replace --lef ./test/ispd18_test1.input.lef --def ./test/ispd18_test1.input.def 

