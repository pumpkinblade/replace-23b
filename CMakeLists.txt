cmake_minimum_required(VERSION 3.12)

project(replace)

set(CMAKE_CXX_STANDARD 17)
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
  ${PROJECT_SOURCE_DIR}/module/CImg_library
)
set(STB_HOME
    ${PROJECT_SOURCE_DIR}/module/stb
)
set(KAHYPAR_HOME
    ${PROJECT_SOURCE_DIR}/module/kahypar
)

add_subdirectory(${DEFLIB_HOME})
add_subdirectory(${LEFLIB_HOME})

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
  src/dreamLegalizer.h

  src/dreamplace/abacus_legalize_cpu.h
  src/dreamplace/bin_assignment_cpu.h
  src/dreamplace/blank.h
  src/dreamplace/compare_cpu.h
  src/dreamplace/function_cpu.h
  src/dreamplace/legalization_db.h
  src/dreamplace/math.h
  src/dreamplace/merge_bin_cpu.h
  src/dreamplace/msg.h
  src/dreamplace/namespace.h
  src/dreamplace/status_summary_cpu.h
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
  src/23bParser.cpp
  src/abacusLegalizer.cpp
  src/partitioner.cpp
  src/macroLegalizer.cpp
  src/terminalModifier.cpp
  src/outputWriter.cpp
  src/abaxLegalizer.cpp
  src/dreamLegalizer.cpp

  src/dreamplace/msg.cpp
  src/dreamplace/greedy_legalize_cpu.cpp
  src/dreamplace/legalize_bin_cpu.cpp
  src/dreamplace/merge_bin_cpu.cpp
)

add_executable(replace ${PLACE_SRC} ${PLACE_HEADER})

add_library( kahypar SHARED IMPORTED )
set_target_properties( kahypar PROPERTIES IMPORTED_LOCATION ${KAHYPAR_HOME}/libkahypar.so )

target_include_directories(replace PRIVATE ${SPDLOG_HOME}/include)
target_include_directories(replace PRIVATE ${TCLAP_HOME}/include)
target_include_directories(replace PRIVATE ${EIGEN_HOME}/include)
target_include_directories(replace PRIVATE ${STB_HOME})
target_include_directories(replace PRIVATE ${CIMG_HOME})

target_link_libraries(replace PRIVATE def)
target_link_libraries(replace PRIVATE lef)
target_link_libraries(replace PRIVATE kahypar)

target_compile_definitions(replace PRIVATE -Dcimg_display=0)
