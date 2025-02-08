cmake_minimum_required(VERSION 3.14)
project(Hist-Tree LANGUAGES CXX)

# C++17 requirement
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

include(CheckCXXSourceRuns)

# Test code checking if AVX2 is available
set(AVX2_TEST_SOURCE "
#include <immintrin.h>
int main() {
    __m256i a = _mm256_set1_epi32(1);
    __m256i b = _mm256_set1_epi32(2);
    __m256i c = _mm256_add_epi32(a, b);
    return _mm256_extract_epi32(c, 0) == 3 ? 0 : 1;
}
")

# temporary AVX2 flag
set(CMAKE_REQUIRED_FLAGS "-mavx2")
check_cxx_source_runs("${AVX2_TEST_SOURCE}" AVX2_SUPPORTED)
set(CMAKE_REQUIRED_FLAGS "")

# Set AVX2 flag if test succeeded
if(AVX2_SUPPORTED)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -mavx2")
    add_definitions(-DUSE_AVX2)
else()
    add_definitions(-DNO_AVX2)
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