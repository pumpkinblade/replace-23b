cmake_minimum_required(VERSION 3.12)

project(replace)

set(CMAKE_CXX_STANDARD 17)
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
set(GLPK_HOME
  ${PROJECT_SOURCE_DIR}/module/glpk
)
set(KAHYPAR_HOME
    ${PROJECT_SOURCE_DIR}/module/kahypar
)
set(MTKAHYPAR_HOME
    ${PROJECT_SOURCE_DIR}/module/mtkahypar
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
  src/plot.cpp
  src/technology.cpp
  src/23bParser.cpp
  src/abacusLegalizer.cpp
  src/partitioner.cpp
  src/macroLegalizer.cpp
  src/outputWriter.cpp
  src/abaxLegalizer.cpp
)

add_executable(replace ${PLACE_SRC} ${PLACE_HEADER})

target_include_directories(replace PRIVATE ${SPDLOG_HOME}/include)
target_include_directories(replace PRIVATE ${TCLAP_HOME}/include)
target_include_directories(replace PRIVATE ${EIGEN_HOME}/include)
target_include_directories(replace PRIVATE ${STB_HOME})
target_include_directories(replace PRIVATE ${CIMG_HOME})
target_include_directories(replace PRIVATE ${GLPK_HOME}/include)
target_include_directories(replace PRIVATE ${KAHYPAR_HOME}/include)

target_compile_definitions(replace PRIVATE -Dcimg_display=0)

if(WIN32)
  target_link_libraries(replace PRIVATE ${GLPK_HOME}/glpk_4_65.lib)
  target_link_libraries(replace PRIVATE ${KAHYPAR_HOME}/kahypar.lib)
  add_custom_command(
    TARGET replace POST_BUILD
    COMMAND ${CMAKE_COMMAND} -E copy 
      ${GLPK_HOME}/glpk_4_65.dll 
      $<TARGET_FILE_DIR:${CMAKE_PROJECT_NAME}>/glpk_4_65.dll
    COMMAND ${CMAKE_COMMAND} -E copy 
      ${KAHYPAR_HOME}/kahypar.dll 
      $<TARGET_FILE_DIR:${CMAKE_PROJECT_NAME}>/kahypar.dll
  )
else()
  add_library( glpk SHARED IMPORTED )
  set_target_properties( glpk PROPERTIES IMPORTED_LOCATION ${GLPK_HOME}/libglpk.so.40 )

  add_library( kahypar SHARED IMPORTED )
  set_target_properties( kahypar PROPERTIES IMPORTED_LOCATION ${KAHYPAR_HOME}/libkahypar.so )

  add_library( mtkahypar SHARED IMPORTED )
  set_target_properties( mtkahypar PROPERTIES IMPORTED_LOCATION ${MTKAHYPAR_HOME}/libmtkahypar.so )

  target_link_libraries(replace PRIVATE kahypar)
  target_link_libraries(replace PRIVATE mtkahypar)

  target_include_directories(replace PRIVATE ${MTKAHYPAR_HOME}/include)
endif()