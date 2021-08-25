cmake_minimum_required(VERSION 3.10)
project(heatmap)

set(CMAKE_CXX_STANDARD 14)

add_executable(heatmap main.cpp gen_data.cpp gen_data.h baseline.cpp baseline.h cell.cpp cell.h CSA.cpp CSA.h utils.cpp utils.h)