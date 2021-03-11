cmake_minimum_required(VERSION 3.17)
project(TEST_LAPACK)

set(CMAKE_CXX_STANDARD 11)

# Flags of sources, headers and libraries
include_directories(/opt/OpenBLAS/include)
link_directories(/opt/OpenBLAS/lib)
set(SOURCE_FILES main.cpp src/Matrix.h src/Matrix.cpp src/Helper.h src/Helper.cpp)

# executable and linking
add_executable(TEST_LAPACK ${SOURCE_FILES})
target_link_libraries(TEST_LAPACK libopenblas.so)