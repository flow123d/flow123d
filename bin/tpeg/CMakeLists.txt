# Project name and languages
project(TPeg CXX)

cmake_minimum_required(VERSION 2.6)

#if (UNIX)
#  message(STATUS "This is UNIX (including OS X and CygWin)")
#endif (UNIX)

set(CMAKE_MODULE_PATH "${CMAKE_ROOT}/Modules")
set(CMAKE_MODULE_PATH "${TPeg_SOURCE_DIR}/CMakeModules" ${CMAKE_MODULE_PATH})

message(STATUS "CMake Module Path: ${CMAKE_MODULE_PATH}")

find_package(Yaml REQUIRED)

add_subdirectory(src)

