# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.10

# Delete rule output on recipe failure.
.DELETE_ON_ERROR:


#=============================================================================
# Special targets provided by cmake.

# Disable implicit rules so canonical targets will work.
.SUFFIXES:


# Remove some rules from gmake that .SUFFIXES does not remove.
SUFFIXES =

.SUFFIXES: .hpux_make_needs_suffix_list


# Suppress display of executed commands.
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
RM = /usr/bin/cmake -E remove -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /home/mti/github-tomerten/IBSlib

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /home/mti/github-tomerten/IBSlib

# Include any dependencies generated for this target.
include CMakeFiles/IBSLib.dir/depend.make

# Include the progress variables for this target.
include CMakeFiles/IBSLib.dir/progress.make

# Include the compile flags for this target's objects.
include CMakeFiles/IBSLib.dir/flags.make

CMakeFiles/IBSLib.dir/libtest.cpp.o: CMakeFiles/IBSLib.dir/flags.make
CMakeFiles/IBSLib.dir/libtest.cpp.o: libtest.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/mti/github-tomerten/IBSlib/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object CMakeFiles/IBSLib.dir/libtest.cpp.o"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/IBSLib.dir/libtest.cpp.o -c /home/mti/github-tomerten/IBSlib/libtest.cpp

CMakeFiles/IBSLib.dir/libtest.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/IBSLib.dir/libtest.cpp.i"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/mti/github-tomerten/IBSlib/libtest.cpp > CMakeFiles/IBSLib.dir/libtest.cpp.i

CMakeFiles/IBSLib.dir/libtest.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/IBSLib.dir/libtest.cpp.s"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/mti/github-tomerten/IBSlib/libtest.cpp -o CMakeFiles/IBSLib.dir/libtest.cpp.s

CMakeFiles/IBSLib.dir/libtest.cpp.o.requires:

.PHONY : CMakeFiles/IBSLib.dir/libtest.cpp.o.requires

CMakeFiles/IBSLib.dir/libtest.cpp.o.provides: CMakeFiles/IBSLib.dir/libtest.cpp.o.requires
	$(MAKE) -f CMakeFiles/IBSLib.dir/build.make CMakeFiles/IBSLib.dir/libtest.cpp.o.provides.build
.PHONY : CMakeFiles/IBSLib.dir/libtest.cpp.o.provides

CMakeFiles/IBSLib.dir/libtest.cpp.o.provides.build: CMakeFiles/IBSLib.dir/libtest.cpp.o


# Object files for target IBSLib
IBSLib_OBJECTS = \
"CMakeFiles/IBSLib.dir/libtest.cpp.o"

# External object files for target IBSLib
IBSLib_EXTERNAL_OBJECTS =

IBSLib: CMakeFiles/IBSLib.dir/libtest.cpp.o
IBSLib: CMakeFiles/IBSLib.dir/build.make
IBSLib: libclog.so
IBSLib: libnumfunc.so
IBSLib: CMakeFiles/IBSLib.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/home/mti/github-tomerten/IBSlib/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Linking CXX executable IBSLib"
	$(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/IBSLib.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
CMakeFiles/IBSLib.dir/build: IBSLib

.PHONY : CMakeFiles/IBSLib.dir/build

CMakeFiles/IBSLib.dir/requires: CMakeFiles/IBSLib.dir/libtest.cpp.o.requires

.PHONY : CMakeFiles/IBSLib.dir/requires

CMakeFiles/IBSLib.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles/IBSLib.dir/cmake_clean.cmake
.PHONY : CMakeFiles/IBSLib.dir/clean

CMakeFiles/IBSLib.dir/depend:
	cd /home/mti/github-tomerten/IBSlib && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/mti/github-tomerten/IBSlib /home/mti/github-tomerten/IBSlib /home/mti/github-tomerten/IBSlib /home/mti/github-tomerten/IBSlib /home/mti/github-tomerten/IBSlib/CMakeFiles/IBSLib.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : CMakeFiles/IBSLib.dir/depend

