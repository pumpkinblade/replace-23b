# $1 : mode : lefdef or 23b
# $2 : case id : int

# set default value
mode=${1:-23b}
testcase=${2:-1}

# clean dumped core file
rm -rf core.[0-9]*

# BUILD
cmake -DCMAKE_BUILD_TYPE=Debug -B build
cmake --build build -j 3
# on build failure, exit
if [ $? -ne 0 ]; then
    echo build failed, exit
    exit 1
fi

# RUN
echo mode: $mode 
rm -rf plot/*/*.jpg
case $mode in 
lefdef)
    ./build/replace --mode ${mode}  --lef ./test/ispd18_test1.input.lef \
                                    --def ./test/ispd18_test1.input.def
    ;;
23b)
    ./build/replace --mode ${mode} --txt23b test/ProblemB_case${testcase}.txt 
    ;;
*)
    echo unknown mode $mode
    exit 1 ;;
esac

