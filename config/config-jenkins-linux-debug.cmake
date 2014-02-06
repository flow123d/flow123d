# Configuration for CI server

set(FLOW_BUILD_TYPE debug)

set(CMAKE_C_COMPILER "/opt/rh/devtoolset-2/root/usr/bin/gcc")
set(CMAKE_CXX_COMPILER "/opt/rh/devtoolset-2/root/usr/bin/g++")

set(CMAKE_VERBOSE_MAKEFILE on)

set(PETSC_INSTALL_CONFIG "bddcml")
set(BDDCML_ROOT "_INSTALL_")
set(USE_PYTHON "yes")

#set(EXTERNAL_PROJECT_DIR "$ENV{HOME}/external_projects/JB_1.7_inputs") 
