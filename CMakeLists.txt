cmake_minimum_required(VERSION 3.27)
project(Sorting_Algorithms)

set(CMAKE_CXX_STANDARD 17)

add_executable(Sorting_Algorithms main.cpp)

find_package(GTest REQUIRED)
include_directories(${GTEST_INCLUDE_DIRS})
target_link_libraries(Sorting_Algorithms ${GTEST_LIBRARIES})
