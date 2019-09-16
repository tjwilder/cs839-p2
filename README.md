# Assignment 1 for CS 839 in F19

## Information
This was built with a specialized build for x64 Windows and build files and instructions are specialized for my machine.

## Building and running
### Cmake to build the solution
cmake . -DPYTHON_LIBRARY:FILEPATH=C:/Python27/libs/python27.lib -DPYTHON_INCLUDE_PATH:FILEPATH=C:/Python27/include

### In VS2017: What you need to do per-project (this could probably be done with cmake, this works):
1. Switch to x64. In Properties > Linker > Command Line, make sure to remove `/machine:X86`
2. In Properties > Linker > General > Additional Library Directories, add K:\Program Files\USD\include\boost-1_65_1\lib64-msvc-14.1
