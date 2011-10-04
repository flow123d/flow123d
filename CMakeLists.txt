# Project name and language(s)
project(Flow123d CXX C)

# CMake version required for pre-build this project
cmake_minimum_required(VERSION 2.6.4)

set(VERSION 1.6.5)
set(CMAKE_BUILD_TYPE release)

# Set all module paths
set(CMAKE_MODULE_PATH "${CMAKE_ROOT}/Modules")
list(APPEND CMAKE_MODULE_PATH "${Flow123d_SOURCE_DIR}/CMake/Modules")
message(STATUS "CMake module path: ${CMAKE_MODULE_PATH}")

set(CMAKE_INCLUDE_CURRENT_DIR ON)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/CMake/makefile.cmake.in
  ${CMAKE_CURRENT_BINARY_DIR}/makefile.in
)
include(${CMAKE_CURRENT_BINARY_DIR}/makefile.in OPTIONAL)

find_package(PETSc REQUIRED)
find_package(PARMETIS REQUIRED)
find_package(Boost REQUIRED)
find_package(Subversion REQUIRED)
find_package(Armadillo REQUIRED)
find_package(GTest)
find_package(CppUnit)

if(Subversion_FOUND)
  Subversion_WC_INFO(${CMAKE_CURRENT_SOURCE_DIR} Flow123d)
  message(STATUS "Current revision is ${Flow123d_WC_REVISION}")
  configure_file(${CMAKE_CURRENT_SOURCE_DIR}/CMake/rev_num.h.in ${CMAKE_CURRENT_SOURCE_DIR}/include/rev_num.h)
endif(Subversion_FOUND)

set(Flow123dIncludeDirectories
    ${CMAKE_SOURCE_DIR}/include
    ${CMAKE_SOURCE_DIR}/include/mesh
    ${CMAKE_SOURCE_DIR}/include/transport
    ${CMAKE_SOURCE_DIR}/src
    ${CMAKE_SOURCE_DIR}/src/flow
    ${CMAKE_SOURCE_DIR}/src/io
    ${CMAKE_SOURCE_DIR}/src/mesh
    ${CMAKE_SOURCE_DIR}/src/reaction
    ${CMAKE_SOURCE_DIR}/src/system
    ${CMAKE_SOURCE_DIR}/src/semchem
    ${CMAKE_SOURCE_DIR}/src/transport
    ${CMAKE_SOURCE_DIR}/third_party/metis-4.0/Lib
    ${PETSC_INCLUDES}
    ${Boost_INCLUDE_DIR}
)

include_directories(${Flow123dIncludeDirectories})

set(Flow123dLinkDirectories
    ${Flow123d_BINARY_dir}/lib
    ${Boost_LIBRARY_DIRS}
)

##############
# set appropriate compiler flags for debug/release compilation modes
#add_definitions("-Wall")

if(CMAKE_BUILD_TYPE STREQUAL "debug")
  # debug mode
  add_definitions("-g -O0")
  add_definitions(-DDEBUG)
else(CMAKE_BUILD_TYPE STREQUAL "debug")
  # Release mode
  add_definitions("-O3")
  add_definitions(-DNDEBUG)
endif(CMAKE_BUILD_TYPE STREQUAL "debug")
##############

add_subdirectory(src)
add_subdirectory(doc)

add_custom_target(benchmark-test 
    COMMAND make testbase
    WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}/tests
)

enable_testing()
add_subdirectory(test_units)

include(CTest)
