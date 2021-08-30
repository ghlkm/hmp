cmake_minimum_required(VERSION 3.10)
project(heatmap)

set(CMAKE_CXX_STANDARD 14)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wl,--no-as-needed -O3 -Wall -ldl")

add_executable(heatmap main.cpp cell.cpp cell.h hmp.cpp hmp.h utils.cpp utils.h vector_operator.h)