cmake -B build
cmake --build build
./build/Debug/replace.exe --lef ./test/ispd18_test1.input.lef --def ./test/ispd18_test1.input.def
./build/Debug/replace.exe --mode 23b --txt23b test/ProblemB_case2.txt 
