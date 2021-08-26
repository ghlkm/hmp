cmake_minimum_required(VERSION 3.10)
project(heatmap)

set(CMAKE_CXX_STANDARD 14)

add_executable(heatmap main.cpp cell.cpp cell.h hmp.cpp hmp.h utils.cpp utils.h vector_operator.h)