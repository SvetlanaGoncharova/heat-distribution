cmake_minimum_required(VERSION 3.6)
project(heat_distribution_lenf)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")

set(SOURCE_FILES
    src/main.cpp
    ../../../initial/INPUT.txt)

add_executable(heat_distribution_lenf ${SOURCE_FILES})