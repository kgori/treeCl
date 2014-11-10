cmake_minimum_required(VERSION 2.8.4)
project(tree_collection)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")

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

add_executable(tree_collection ${SOURCE_FILES} src/TreeCollection.cc)
add_executable(test ${SOURCE_FILES} test.cpp)
