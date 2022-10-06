cmake_minimum_required(VERSION 3.23)
project(SchroedingerSolverTemplate)

set(CMAKE_CXX_STANDARD 20)

find_package(Boost)
include_directories(${Boost_INCLUDE_DIR})

add_executable(SchroedingerSolverTemplate main.cpp)
