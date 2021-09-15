cmake_minimum_required(VERSION 3.10)
project(heatmap)

set(CMAKE_CXX_STANDARD 14)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 -Wall -ldl")

include_directories(${CMAKE_SOURCE_DIR}/rtree)

add_executable(heatmap main.cpp cell.cpp cell.h hmp.cpp hmp.h utils.cpp utils.h vector_operator.h
        rtree/collection.h
        rtree/collection.cpp
        rtree/filemem.h
        rtree/filemem.cpp
        rtree/global.h
        rtree/global.cpp
        rtree/header.h
        rtree/hypercube.h
        rtree/hypercube.cpp
        rtree/mem.h
        rtree/param.h
        rtree/param.cpp
        rtree/point.h
        rtree/point.cpp
        rtree/rentry.h
        rtree/rentry.cpp
        rtree/rnode.h
        rtree/rnode.cpp
        rtree/rtree.h
        rtree/rtree.cpp
        rtree/tgs.h
        rtree/tgs.cpp
        rtree/virtualRNode.h
        rtree/virtualRNode.cpp)