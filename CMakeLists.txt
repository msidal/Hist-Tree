# minimum version of CMake
cmake_minimum_required(VERSION 3.14)

# projectname and language
project(Hist-Tree LANGUAGES CXX)

# C++-standard
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# GoogleTest integration via FetchContent
#TODO

# include directories
include_directories(include)

# source files
file(GLOB SOURCES src/*.cpp)

# define executable
add_executable(hist-tree main.cpp ${SOURCES})

