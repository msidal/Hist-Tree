# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.28

# Delete rule output on recipe failure.
.DELETE_ON_ERROR:

#=============================================================================
# Special targets provided by cmake.

# Disable implicit rules so canonical targets will work.
.SUFFIXES:

# Disable VCS-based implicit rules.
% : %,v

# Disable VCS-based implicit rules.
% : RCS/%

# Disable VCS-based implicit rules.
% : RCS/%,v

# Disable VCS-based implicit rules.
% : SCCS/s.%

# Disable VCS-based implicit rules.
% : s.%

.SUFFIXES: .hpux_make_needs_suffix_list

# Command-line flag to silence nested $(MAKE).
$(VERBOSE)MAKESILENT = -s

#Suppress display of executed commands.
$(VERBOSE).SILENT:

# A target that is always out of date.
cmake_force:
.PHONY : cmake_force

#=============================================================================
# Set environment variables for the build.

# The shell in which to execute make rules.
SHELL = /bin/sh

# The CMake executable.
CMAKE_COMMAND = /usr/bin/cmake

# The command to remove a file.
RM = /usr/bin/cmake -E rm -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /home/mert/Hist-Tree

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /home/mert/Hist-Tree/build

# Include any dependencies generated for this target.
include CMakeFiles/hist-tree.dir/depend.make
# Include any dependencies generated by the compiler for this target.
include CMakeFiles/hist-tree.dir/compiler_depend.make

# Include the progress variables for this target.
include CMakeFiles/hist-tree.dir/progress.make

# Include the compile flags for this target's objects.
include CMakeFiles/hist-tree.dir/flags.make

CMakeFiles/hist-tree.dir/main.cpp.o: CMakeFiles/hist-tree.dir/flags.make
CMakeFiles/hist-tree.dir/main.cpp.o: /home/mert/Hist-Tree/main.cpp
CMakeFiles/hist-tree.dir/main.cpp.o: CMakeFiles/hist-tree.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green --progress-dir=/home/mert/Hist-Tree/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object CMakeFiles/hist-tree.dir/main.cpp.o"
	/usr/bin/clang++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -MD -MT CMakeFiles/hist-tree.dir/main.cpp.o -MF CMakeFiles/hist-tree.dir/main.cpp.o.d -o CMakeFiles/hist-tree.dir/main.cpp.o -c /home/mert/Hist-Tree/main.cpp

CMakeFiles/hist-tree.dir/main.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green "Preprocessing CXX source to CMakeFiles/hist-tree.dir/main.cpp.i"
	/usr/bin/clang++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/mert/Hist-Tree/main.cpp > CMakeFiles/hist-tree.dir/main.cpp.i

CMakeFiles/hist-tree.dir/main.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green "Compiling CXX source to assembly CMakeFiles/hist-tree.dir/main.cpp.s"
	/usr/bin/clang++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/mert/Hist-Tree/main.cpp -o CMakeFiles/hist-tree.dir/main.cpp.s

# Object files for target hist-tree
hist__tree_OBJECTS = \
"CMakeFiles/hist-tree.dir/main.cpp.o"

# External object files for target hist-tree
hist__tree_EXTERNAL_OBJECTS =

hist-tree: CMakeFiles/hist-tree.dir/main.cpp.o
hist-tree: CMakeFiles/hist-tree.dir/build.make
hist-tree: CMakeFiles/hist-tree.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green --bold --progress-dir=/home/mert/Hist-Tree/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Linking CXX executable hist-tree"
	$(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/hist-tree.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
CMakeFiles/hist-tree.dir/build: hist-tree
.PHONY : CMakeFiles/hist-tree.dir/build

CMakeFiles/hist-tree.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles/hist-tree.dir/cmake_clean.cmake
.PHONY : CMakeFiles/hist-tree.dir/clean

CMakeFiles/hist-tree.dir/depend:
	cd /home/mert/Hist-Tree/build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/mert/Hist-Tree /home/mert/Hist-Tree /home/mert/Hist-Tree/build /home/mert/Hist-Tree/build /home/mert/Hist-Tree/build/CMakeFiles/hist-tree.dir/DependInfo.cmake "--color=$(COLOR)"
.PHONY : CMakeFiles/hist-tree.dir/depend

