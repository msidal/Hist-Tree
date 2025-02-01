# minimum version of CMake
cmake_minimum_required(VERSION 3.14)

# projectname and language
project(Hist-Tree LANGUAGES CXX)

# C++-standard
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)


#flags for optimized build
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -mavx2 -mfma -O3 -Wall -Wextra -Wpedantic -Werror -fopenmp")

# Debugging flags
#set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -mavx2 -mfma -Og -ggdb -fno-omit-frame-pointer")

# GoogleTest and GoogleBenchmark integration via FetchContent
include(FetchContent)

FetchContent_Declare(
  googletest
  URL https://github.com/google/googletest/archive/03597a01ee50ed33e9dfd640b249b4be3799d395.zip
)

FetchContent_Declare(
  googlebenchmark
  URL https://github.com/google/benchmark/archive/refs/tags/v1.8.2.zip
)

FetchContent_MakeAvailable(googletest googlebenchmark)

# include directories
include_directories(include)

# source files if necessary
# file(GLOB SOURCES src/*.cpp)

# define executable
add_executable(hist-tree main.cpp) # ${SOURCES})

# Tests
enable_testing()

add_executable(
  builder_test
  tests/builder_test.cpp
)

target_link_libraries(
  builder_test
  GTest::gtest_main
)

add_executable(
  histtree_test
  tests/histtree_test.cpp
)

target_link_libraries(
  histtree_test
  GTest::gtest_main
)

include(GoogleTest)
gtest_discover_tests(builder_test)
gtest_discover_tests(histtree_test)

# Benchmarks
add_executable(
  histtree_benchmark
  benchmarks/benchmarkTest.cpp
)

target_link_libraries(
  histtree_benchmark
  benchmark::benchmark
)
