cmake_minimum_required(VERSION 3.17)
project(CPPKIT)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fconcepts")

# Flags of sources, headers and libraries
set_directory_properties(PROPERTIES COMPILE_DEFINITIONS_DEBUG "_DEBUG")
include_directories(/opt/OpenBLAS/include)
link_directories(/opt/OpenBLAS/lib)
set(SOURCE_FILES main.cpp src/Matrix.h src/Matrix.cpp src/Helper.h src/Helper.cpp src/LAPACK.h src/BLAS.h src/Vector.h src/Vector.cpp)

# executable and linking
add_executable(CPPKIT ${SOURCE_FILES})
target_link_libraries(CPPKIT libopenblas.so)