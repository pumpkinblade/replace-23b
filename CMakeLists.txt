cmake_minimum_required(VERSION 3.12)

project(replace)

# CIMG is a header only package
option(USE_CIMG "use cimg" OFF)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

set(REPLACE_HOME ${PROJECT_SOURCE_DIR} )
set(DEFLIB_HOME
    ${PROJECT_SOURCE_DIR}/module/def
)
set(LEFLIB_HOME
    ${PROJECT_SOURCE_DIR}/module/lef
)
set(SPDLOG_HOME
    ${PROJECT_SOURCE_DIR}/module/spdlog
)
set(TCLAP_HOME
    ${PROJECT_SOURCE_DIR}/module/tclap
)
set(EIGEN_HOME
    ${PROJECT_SOURCE_DIR}/module/eigen
)
set(CIMG_HOME
    ${PROJECT_SOURCE_DIR}/module/CImg
)

add_subdirectory(${DEFLIB_HOME})
add_subdirectory(${LEFLIB_HOME})

if(USE_CIMG)
    find_package(JPEG REQUIRED)
    find_package(Threads REQUIRED)
    #find_package(Cimg CONFIG REQUIRED)
endif()

set(PLACE_HEADER
  src/log.h
  src/parser.h
  src/placerBase.h
  src/replace.h
  src/initialPlace.h
  src/fft.h
  src/nesterovBase.h
  src/nesterovPlace.h
  src/point.h
  src/plot.h
  src/technology.h
  src/placer23b.h
  src/abacusLegalizer.h
)

set(PLACE_SRC
  src/main.cpp
  src/log.cpp
  src/lefdefParser.cpp
  src/placerBase.cpp
  src/replace.cpp
  src/initialPlace.cpp
  src/fft.cpp
  src/fftsg.cpp
  src/fftsg2d.cpp
  src/fftsg3d.cpp
  src/nesterovBase.cpp
  src/nesterovPlace.cpp
  src/point.cpp
  src/plot.cpp
  src/technology.cpp
  src/placer23b.cpp
  src/23bParser.cpp
  src/abacusLegalizer.cpp
)

add_executable(replace ${PLACE_SRC} ${PLACE_HEADER})

target_include_directories(replace PRIVATE ${SPDLOG_HOME}/include)
target_include_directories(replace PRIVATE ${TCLAP_HOME}/include)
target_include_directories(replace PRIVATE ${EIGEN_HOME}/include)
if(USE_CIMG)
    target_include_directories(replace PRIVATE ${CIMG_HOME})
endif()

target_link_libraries(replace PRIVATE def)
target_link_libraries(replace PRIVATE lef)

if(USE_CIMG)
    # target_link_libraries(replace PRIVATE CImg::CImg)
    target_link_libraries(replace PRIVATE ${JPEG_LIBRARIES})
    target_link_libraries(replace PRIVATE Threads::Threads)
    target_compile_options(replace PRIVATE -DENABLE_CIMG_LIB)
    target_compile_options(replace PRIVATE -Dcimg_use_jpeg)
endif()
