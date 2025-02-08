cmake_minimum_required(VERSION 3.14)
project(Hist-Tree LANGUAGES CXX)

# C++17 requirement
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Check for AVX2 support
include(CheckCXXCompilerFlag)
check_cxx_compiler_flag("-mavx2" COMPILER_SUPPORTS_AVX2)
if(COMPILER_SUPPORTS_AVX2)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -mavx2 -mfma")
endif()

# Compiler flags
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 -Wall -Wextra -Wpedantic -Werror")

# OpenMP support
find_package(OpenMP)
if(OpenMP_CXX_FOUND)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
endif()

# Configure FetchContent
include(FetchContent)

# Configure GoogleTest
FetchContent_Declare(
    googletest
    URL https://github.com/google/googletest/archive/03597a01ee50ed33e9dfd640b249b4be3799d395.zip
    OVERRIDE_FIND_PACKAGE
)

# Configure Google Benchmark
FetchContent_Declare(
    googlebenchmark
    URL https://github.com/google/benchmark/archive/refs/tags/v1.8.2.zip
    OVERRIDE_FIND_PACKAGE
)

# Make dependencies available
FetchContent_MakeAvailable(googletest googlebenchmark)

# Include directories
include_directories(include)

# Main executable
add_executable(hist-tree main.cpp)

# Testing configuration
enable_testing()

# Test executables
add_executable(builder_test tests/builder_test.cpp)
target_link_libraries(builder_test GTest::gtest_main)

add_executable(histtree_test tests/histtree_test.cpp)
target_link_libraries(histtree_test GTest::gtest_main)

# Discover tests
include(GoogleTest)
gtest_discover_tests(builder_test)
gtest_discover_tests(histtree_test)

# Benchmark executable
add_executable(bench benchmarks/benchmark.cpp)
target_link_libraries(bench benchmark::benchmark)