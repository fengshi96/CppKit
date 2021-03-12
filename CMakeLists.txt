cmake_minimum_required(VERSION 3.17)
project(CPPKIT)

set(CMAKE_CXX_STANDARD 20)

# Flags of sources, headers and libraries
include_directories(/opt/OpenBLAS/include)
link_directories(/opt/OpenBLAS/lib)
set(SOURCE_FILES main.cpp src/Matrix.h src/Matrix.cpp src/Helper.h src/Helper.cpp src/LAPACK.h src/BLAS.h)

# executable and linking
add_executable(CPPKIT ${SOURCE_FILES})
target_link_libraries(CPPKIT libopenblas.so)