cmake_minimum_required(VERSION 2.8.4)
project(tree_collection)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")
include_directories(src/eigen3)
set(SOURCE_FILES
    cython/wrapper.cpp
    src/MinSqTree.cc
    src/MinSqTree.h
    src/newick.cc
    src/newick.h
    src/PhyTree.h
    src/ProblemParser.cc
    src/ProblemParser.h
    )


add_executable(tst ${SOURCE_FILES} test.cpp)
