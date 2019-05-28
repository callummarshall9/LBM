# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.14

# Default target executed when no arguments are given to make.
default_target: all

.PHONY : default_target

# Allow only one "make -f Makefile2" at a time, but pass parallelism.
.NOTPARALLEL:


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
CMAKE_SOURCE_DIR = /home/callummarshall/Documents/Development/Cpp/LBM

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /home/callummarshall/Documents/Development/Cpp/LBM

#=============================================================================
# Targets provided globally by CMake.

# Special rule for the target rebuild_cache
rebuild_cache:
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --cyan "Running CMake to regenerate build system..."
	/usr/bin/cmake -S$(CMAKE_SOURCE_DIR) -B$(CMAKE_BINARY_DIR)
.PHONY : rebuild_cache

# Special rule for the target rebuild_cache
rebuild_cache/fast: rebuild_cache

.PHONY : rebuild_cache/fast

# Special rule for the target edit_cache
edit_cache:
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --cyan "Running CMake cache editor..."
	/usr/bin/ccmake -S$(CMAKE_SOURCE_DIR) -B$(CMAKE_BINARY_DIR)
.PHONY : edit_cache

# Special rule for the target edit_cache
edit_cache/fast: edit_cache

.PHONY : edit_cache/fast

# The main all target
all: cmake_check_build_system
	$(CMAKE_COMMAND) -E cmake_progress_start /home/callummarshall/Documents/Development/Cpp/LBM/CMakeFiles /home/callummarshall/Documents/Development/Cpp/LBM/CMakeFiles/progress.marks
	$(MAKE) -f CMakeFiles/Makefile2 all
	$(CMAKE_COMMAND) -E cmake_progress_start /home/callummarshall/Documents/Development/Cpp/LBM/CMakeFiles 0
.PHONY : all

# The main clean target
clean:
	$(MAKE) -f CMakeFiles/Makefile2 clean
.PHONY : clean

# The main clean target
clean/fast: clean

.PHONY : clean/fast

# Prepare targets for installation.
preinstall: all
	$(MAKE) -f CMakeFiles/Makefile2 preinstall
.PHONY : preinstall

# Prepare targets for installation.
preinstall/fast:
	$(MAKE) -f CMakeFiles/Makefile2 preinstall
.PHONY : preinstall/fast

# clear depends
depend:
	$(CMAKE_COMMAND) -S$(CMAKE_SOURCE_DIR) -B$(CMAKE_BINARY_DIR) --check-build-system CMakeFiles/Makefile.cmake 1
.PHONY : depend

#=============================================================================
# Target rules for targets named main_cmdline

# Build rule for target.
main_cmdline: cmake_check_build_system
	$(MAKE) -f CMakeFiles/Makefile2 main_cmdline
.PHONY : main_cmdline

# fast build rule for target.
main_cmdline/fast:
	$(MAKE) -f CMakeFiles/main_cmdline.dir/build.make CMakeFiles/main_cmdline.dir/build
.PHONY : main_cmdline/fast

#=============================================================================
# Target rules for targets named main_gui

# Build rule for target.
main_gui: cmake_check_build_system
	$(MAKE) -f CMakeFiles/Makefile2 main_gui
.PHONY : main_gui

# fast build rule for target.
main_gui/fast:
	$(MAKE) -f CMakeFiles/main_gui.dir/build.make CMakeFiles/main_gui.dir/build
.PHONY : main_gui/fast

LBM/source/LBM.o: LBM/source/LBM.cpp.o

.PHONY : LBM/source/LBM.o

# target to build an object file
LBM/source/LBM.cpp.o:
	$(MAKE) -f CMakeFiles/main_cmdline.dir/build.make CMakeFiles/main_cmdline.dir/LBM/source/LBM.cpp.o
	$(MAKE) -f CMakeFiles/main_gui.dir/build.make CMakeFiles/main_gui.dir/LBM/source/LBM.cpp.o
.PHONY : LBM/source/LBM.cpp.o

LBM/source/LBM.i: LBM/source/LBM.cpp.i

.PHONY : LBM/source/LBM.i

# target to preprocess a source file
LBM/source/LBM.cpp.i:
	$(MAKE) -f CMakeFiles/main_cmdline.dir/build.make CMakeFiles/main_cmdline.dir/LBM/source/LBM.cpp.i
	$(MAKE) -f CMakeFiles/main_gui.dir/build.make CMakeFiles/main_gui.dir/LBM/source/LBM.cpp.i
.PHONY : LBM/source/LBM.cpp.i

LBM/source/LBM.s: LBM/source/LBM.cpp.s

.PHONY : LBM/source/LBM.s

# target to generate assembly for a file
LBM/source/LBM.cpp.s:
	$(MAKE) -f CMakeFiles/main_cmdline.dir/build.make CMakeFiles/main_cmdline.dir/LBM/source/LBM.cpp.s
	$(MAKE) -f CMakeFiles/main_gui.dir/build.make CMakeFiles/main_gui.dir/LBM/source/LBM.cpp.s
.PHONY : LBM/source/LBM.cpp.s

LBM/source/vector3.o: LBM/source/vector3.cpp.o

.PHONY : LBM/source/vector3.o

# target to build an object file
LBM/source/vector3.cpp.o:
	$(MAKE) -f CMakeFiles/main_cmdline.dir/build.make CMakeFiles/main_cmdline.dir/LBM/source/vector3.cpp.o
	$(MAKE) -f CMakeFiles/main_gui.dir/build.make CMakeFiles/main_gui.dir/LBM/source/vector3.cpp.o
.PHONY : LBM/source/vector3.cpp.o

LBM/source/vector3.i: LBM/source/vector3.cpp.i

.PHONY : LBM/source/vector3.i

# target to preprocess a source file
LBM/source/vector3.cpp.i:
	$(MAKE) -f CMakeFiles/main_cmdline.dir/build.make CMakeFiles/main_cmdline.dir/LBM/source/vector3.cpp.i
	$(MAKE) -f CMakeFiles/main_gui.dir/build.make CMakeFiles/main_gui.dir/LBM/source/vector3.cpp.i
.PHONY : LBM/source/vector3.cpp.i

LBM/source/vector3.s: LBM/source/vector3.cpp.s

.PHONY : LBM/source/vector3.s

# target to generate assembly for a file
LBM/source/vector3.cpp.s:
	$(MAKE) -f CMakeFiles/main_cmdline.dir/build.make CMakeFiles/main_cmdline.dir/LBM/source/vector3.cpp.s
	$(MAKE) -f CMakeFiles/main_gui.dir/build.make CMakeFiles/main_gui.dir/LBM/source/vector3.cpp.s
.PHONY : LBM/source/vector3.cpp.s

cmdline/main.o: cmdline/main.cpp.o

.PHONY : cmdline/main.o

# target to build an object file
cmdline/main.cpp.o:
	$(MAKE) -f CMakeFiles/main_cmdline.dir/build.make CMakeFiles/main_cmdline.dir/cmdline/main.cpp.o
.PHONY : cmdline/main.cpp.o

cmdline/main.i: cmdline/main.cpp.i

.PHONY : cmdline/main.i

# target to preprocess a source file
cmdline/main.cpp.i:
	$(MAKE) -f CMakeFiles/main_cmdline.dir/build.make CMakeFiles/main_cmdline.dir/cmdline/main.cpp.i
.PHONY : cmdline/main.cpp.i

cmdline/main.s: cmdline/main.cpp.s

.PHONY : cmdline/main.s

# target to generate assembly for a file
cmdline/main.cpp.s:
	$(MAKE) -f CMakeFiles/main_cmdline.dir/build.make CMakeFiles/main_cmdline.dir/cmdline/main.cpp.s
.PHONY : cmdline/main.cpp.s

gui/main.o: gui/main.cpp.o

.PHONY : gui/main.o

# target to build an object file
gui/main.cpp.o:
	$(MAKE) -f CMakeFiles/main_gui.dir/build.make CMakeFiles/main_gui.dir/gui/main.cpp.o
.PHONY : gui/main.cpp.o

gui/main.i: gui/main.cpp.i

.PHONY : gui/main.i

# target to preprocess a source file
gui/main.cpp.i:
	$(MAKE) -f CMakeFiles/main_gui.dir/build.make CMakeFiles/main_gui.dir/gui/main.cpp.i
.PHONY : gui/main.cpp.i

gui/main.s: gui/main.cpp.s

.PHONY : gui/main.s

# target to generate assembly for a file
gui/main.cpp.s:
	$(MAKE) -f CMakeFiles/main_gui.dir/build.make CMakeFiles/main_gui.dir/gui/main.cpp.s
.PHONY : gui/main.cpp.s

gui/source/worker.o: gui/source/worker.cpp.o

.PHONY : gui/source/worker.o

# target to build an object file
gui/source/worker.cpp.o:
	$(MAKE) -f CMakeFiles/main_gui.dir/build.make CMakeFiles/main_gui.dir/gui/source/worker.cpp.o
.PHONY : gui/source/worker.cpp.o

gui/source/worker.i: gui/source/worker.cpp.i

.PHONY : gui/source/worker.i

# target to preprocess a source file
gui/source/worker.cpp.i:
	$(MAKE) -f CMakeFiles/main_gui.dir/build.make CMakeFiles/main_gui.dir/gui/source/worker.cpp.i
.PHONY : gui/source/worker.cpp.i

gui/source/worker.s: gui/source/worker.cpp.s

.PHONY : gui/source/worker.s

# target to generate assembly for a file
gui/source/worker.cpp.s:
	$(MAKE) -f CMakeFiles/main_gui.dir/build.make CMakeFiles/main_gui.dir/gui/source/worker.cpp.s
.PHONY : gui/source/worker.cpp.s

# Help Target
help:
	@echo "The following are some of the valid targets for this Makefile:"
	@echo "... all (the default if no target is provided)"
	@echo "... clean"
	@echo "... depend"
	@echo "... rebuild_cache"
	@echo "... edit_cache"
	@echo "... main_cmdline"
	@echo "... main_gui"
	@echo "... LBM/source/LBM.o"
	@echo "... LBM/source/LBM.i"
	@echo "... LBM/source/LBM.s"
	@echo "... LBM/source/vector3.o"
	@echo "... LBM/source/vector3.i"
	@echo "... LBM/source/vector3.s"
	@echo "... cmdline/main.o"
	@echo "... cmdline/main.i"
	@echo "... cmdline/main.s"
	@echo "... gui/main.o"
	@echo "... gui/main.i"
	@echo "... gui/main.s"
	@echo "... gui/source/worker.o"
	@echo "... gui/source/worker.i"
	@echo "... gui/source/worker.s"
.PHONY : help



#=============================================================================
# Special targets to cleanup operation of make.

# Special rule to run CMake to check the build system integrity.
# No rule that depends on this can have commands that come from listfiles
# because they might be regenerated.
cmake_check_build_system:
	$(CMAKE_COMMAND) -S$(CMAKE_SOURCE_DIR) -B$(CMAKE_BINARY_DIR) --check-build-system CMakeFiles/Makefile.cmake 0
.PHONY : cmake_check_build_system
