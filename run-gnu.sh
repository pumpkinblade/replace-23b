# $1 : case id : int

# set default value
testcase=${1:-1}

# clean dumped core file
rm -rf core.[0-9]*

# BUILD
cmake -DCMAKE_BUILD_TYPE=Debug -B build
cmake --build build -j 4
# on build failure, exit
if [ $? -ne 0 ]; then
    echo build failed, exit
    exit 1
fi

# RUN
rm -rf plot/*/*.jpg
./build/Debug/replace -i test/ProblemB_case${testcase}.txt -o result/ProblemB_case${testcase}_result.txt

