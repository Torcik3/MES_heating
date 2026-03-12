cmake_minimum_required(VERSION 4.1)
project(MES)

set(CMAKE_CXX_STANDARD 20)

add_executable(MES main.cpp
        mesh.cpp
        femCalculator.cpp
        femCalculator.h)
