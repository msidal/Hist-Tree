# minimum version of CMake
cmake_minimum_required(VERSION 3.14)

# projectname and language
project(Hist-Tree LANGUAGES CXX)

# C++-standard
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)


#flags for SIMD
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -mavx2 -mfma")

# GoogleTest integration via FetchContent
include(FetchContent)
FetchContent_Declare(
  googletest
  URL https://github.com/google/googletest/archive/03597a01ee50ed33e9dfd640b249b4be3799d395.zip
)

FetchContent_MakeAvailable(googletest)

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
  test/builder_test.cpp
)

target_link_libraries(
  builder_test
  GTest::gtest_main
)

add_executable(
  histtree_test
  test/histtree_test.cpp
)

target_link_libraries(
  histtree_test
  GTest::gtest_main
)

include(GoogleTest)
gtest_discover_tests(builder_test)
gtest_discover_tests(histtree_test)


