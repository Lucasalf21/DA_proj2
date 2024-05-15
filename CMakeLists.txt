cmake_minimum_required(VERSION 3.28)
project(proj2_2324)

set(CMAKE_CXX_STANDARD 17)

add_executable(proj2_2324 src/main.cpp
        src/data_structures/Graph.h
        src/data_structures/MutablePriorityQueue.h
        src/data_structures/UFDS.cpp
        src/GraphLoader.h
        src/GraphLoader.cpp)
