cmake_minimum_required(VERSION 2.8)

project(arma_iterative)

set(CMAKE_C_COMPILER "/usr/bin/gcc")
set(CMAKE_CXX_COMPILER "/usr/bin/g++")
set(CMAKE_CXX_FLAGS "-O3")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")
set(CMAKE_INCLUDE_CURRENT_DIR ON)

find_package(OpenMP)
find_package(Armadillo REQUIRED)
find_package(BLAS REQUIRED)
find_package(LAPACK REQUIRED)
#find_package(ARPACK REQUIRED)

include_directories("armadillo")
include_directories("armadillo/armadillo_bits")
add_subdirectory(armadillo/)
include_directories("lis")
add_subdirectory(lis/)
include_directories("lis/include")
add_subdirectory(lis/include)

set(INSTALLED_LIBRARY openblas lapack superlu arpack)

add_library(${PROJECT_NAME} iterative_solve.hpp)
if(OPENMP_FOUND)
    target_link_libraries(${PROJECT_NAME} armadillo lis OpenMP::OpenMP_CXX ${INSTALLED_LIBRARY})
else()
    target_link_libraries(${PROJECT_NAME} armadillo lis ${INSTALLED_LIBRARY})
endif()
set_target_properties(${PROJECT_NAME} PROPERTIES LINKER_LANGUAGE CXX)

add_executable(test3 arma_iterative.cpp)
target_link_libraries(test3 arma_iterative)

