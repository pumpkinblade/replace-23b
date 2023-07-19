# $1 : mode : lefdef or 23b, default to lefdef

mode=${1:-lefdef}

cmake -DCMAKE_BUILD_TYPE=Debug -B build
cmake --build build
# on build failure, exit
if [ $? -ne 0 ]; then
    echo build failed, exit
    exit 1
fi

echo mode: $mode 
rm -rf plot/*/*.jpg
case $mode in 
lefdef)
    ./build/replace --mode ${mode}  --lef ./test/ispd18_test1.input.lef \
                                    --def ./test/ispd18_test1.input.def
    ;;
23b)
    ./build/replace --mode ${mode} --txt23b test/ProblemB_case1.txt 
    ;;
*)
    echo unknown mode $mode
    exit 1 ;;
esac

