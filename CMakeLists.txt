cmake_minimum_required(VERSION 3.12)

project(replace)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DUSE_CIMG=ON -lpthread ")
SET(CMAKE_CXX_FLAGS_DEBUG "$ENV{CXXFLAGS} -O0 -Wall -g2 -ggdb -DUSE_CIMG=ON")
SET(CMAKE_CXX_FLAGS_RELEASE "$ENV{CXXFLAGS} -O3 -Wall -DUSE_CIMG=ON")
# option(USE_CIMG "use cimg" ON)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)



set(THREADS_PREFER_PTHREAD_FLAG ON)

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
# add_subdirectory(${CIMG_HOME})

if(USE_CIMG)
    find_package(JPEG REQUIRED)
    # find_package(Cimg CONFIG REQUIRED)
    find_package(X11 REQUIRED)
    add_definitions(-DENABLE_CIMG_LIB=1)
    set (CMAKE_SHARED_LINKER_FLAGS "${CMAKE_SHARED_LINKER_FLAGS} ${X11_LIBRARIES} ${JPEG_LIBRARIES}")
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
  src/partitioner.h
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
  src/partitioner.cpp
)

add_executable(replace ${PLACE_SRC} ${PLACE_HEADER})

target_include_directories(replace PRIVATE ${SPDLOG_HOME}/include)
target_include_directories(replace PRIVATE ${TCLAP_HOME}/include)
target_include_directories(replace PRIVATE ${EIGEN_HOME}/include)

target_link_libraries(replace PRIVATE def)
target_link_libraries(replace PRIVATE lef)

if(USE_CIMG)
  target_link_libraries( replace
    PRIVATE
    ${X11_LIBRARIES} 
    ${JPEG_LIBRARIES} 
  )
  target_include_directories( replace
    PUBLIC
    # For CImg
    ${REPLACE_HOME}/module/CImg
    ${X11_INCLUDE_DIR}
    ${JPEG_INCLUDE_DIR}
  )
endif()

# if(USE_CIMG)
#     target_link_libraries(replace PRIVATE CImg::CImg)
#     target_link_libraries(replace PRIVATE ${JPEG_LIBRARIES})
#     target_compile_options(replace PRIVATE -DENABLE_CIMG_LIB)
#     target_compile_options(replace PRIVATE -Dcimg_use_jpeg)
# endif()
