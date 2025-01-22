#!/bin/bash

# Simple build and test script for hist-tree
# you might have to run chmod +x test.sh to make it executable

echo "Creating the build directory if it doesn't exist..."
mkdir -p build

echo "Switching to the build directory..."
cd build || exit 1

echo "Running CMake..."
cmake .. || { echo "CMake failed!"; exit 1; }

echo "Building the project with Make..."
make || { echo "Make failed!"; exit 1; }

echo "Running the tests (ctest)..."
ctest|| { echo "Test execution failed!"; exit 1; }