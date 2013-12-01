# Configuration for CI server

set(CMAKE_C_COMPILER "/usr/bin/gcc")
set(CMAKE_CXX_COMPILER "/usr/bin/g++")

set(FLOW_CC_FLAGS "-O3 -DDEBUG_PROFILER")


set(CMAKE_VERBOSE_MAKEFILE on)

set(INSTALL_PETSC_BDDCML "yes")
set(USE_PYTHON "yes")
