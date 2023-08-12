cmake_minimum_required(VERSION 3.12)

project(replace)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

set(REPLACE_HOME ${PROJECT_SOURCE_DIR} )
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
  ${PROJECT_SOURCE_DIR}/module/CImg_library
)
set(STB_HOME
  ${PROJECT_SOURCE_DIR}/module/stb
)

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
  src/abacusLegalizer.h
  src/partitioner.h
  src/macroLegalizer.h
  src/terminalModifier.h
  src/outputWriter.h
  src/abaxLegalizer.h
)

set(PLACE_SRC
  src/main.cpp
  src/log.cpp
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
  src/23bParser.cpp
  src/abacusLegalizer.cpp
  src/partitioner.cpp
  src/macroLegalizer.cpp
  src/terminalModifier.cpp
  src/outputWriter.cpp
  src/abaxLegalizer.cpp
)

add_executable(replace ${PLACE_SRC} ${PLACE_HEADER})

target_include_directories(replace PRIVATE ${SPDLOG_HOME}/include)
target_include_directories(replace PRIVATE ${TCLAP_HOME}/include)
target_include_directories(replace PRIVATE ${EIGEN_HOME}/include)
target_include_directories(replace PRIVATE ${STB_HOME})
target_include_directories(replace PRIVATE ${CIMG_HOME})

target_compile_definitions(replace PRIVATE -Dcimg_display=0)
