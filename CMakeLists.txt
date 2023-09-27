cmake_minimum_required(VERSION 3.10)

project(concepts VERSION 1.0)

# find_package (Eigen3 REQUIRED NO_MODULE)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED True)

set(hdrs 
inc/point_set.hpp
inc/operators.hpp
inc/groups.hpp
inc/fields.hpp
)

set(srcs
src/main.cc
)

add_executable(concepts ${hdrs} ${srcs})
target_include_directories(concepts PUBLIC inc eigen-3.4.0)