# CMAKE system TODO
#
#
# * report compiler flags and used libraries in a summarizing table before own build
#
# * use CTest to perform unit tests
# 
# * use CPack to  make Windows and Linux packages
#
#

# check build_tree link,
# this do not duplicate of the similar call in main makefile,
# since CMakeList.txt can be called by generated makefiles.
execute_process(COMMAND bin/git_post_checkout_hook WORKING_DIRECTORY ${CMAKE_SOURCE_DIR})



# Include settings given by user first
set(CONFIG_CMAKE "config.cmake")
set(CONFIG_CMAKE_DEFAULT "config.default.cmake")

# require default config which is a part of a repo
MESSAGE(STATUS "Including config file: ${CMAKE_SOURCE_DIR}/${CONFIG_CMAKE_DEFAULT}")
include("${CMAKE_SOURCE_DIR}/${CONFIG_CMAKE_DEFAULT}")

# try to read config.cmake if exists
if (EXISTS "${CMAKE_SOURCE_DIR}/${CONFIG_CMAKE}")
    MESSAGE(STATUS "Including config file: ${CMAKE_SOURCE_DIR}/${CONFIG_CMAKE}")
    include("${CMAKE_SOURCE_DIR}/${CONFIG_CMAKE}")
endif()

# if FLOW_BUILD_TYPE is AUTO, it means we did not loaded config.cmake file
# if this is running it our docker container, it'll be fine, otherwise we fail
if (FLOW_BUILD_TYPE STREQUAL "AUTO")
  SET(DOCKER_VERSION_FILE "/.dockerversion")

  # read dockerversion file
  if (EXISTS "${DOCKER_VERSION_FILE}")
    FILE(READ "${DOCKER_VERSION_FILE}" DOCKER_VERSION)
    string(STRIP "${DOCKER_VERSION}" DOCKER_VERSION)
    MESSAGE(STATUS "DOCKER_VERSION=${DOCKER_VERSION}")
  else()
    MESSAGE(FATAL_ERROR "Not running in docker, need configuration file: ${CONFIG_CMAKE}. See config.cmake.template.")
  endif()
  
  # try to detect debug/release
  if ("${DOCKER_VERSION}" MATCHES "dbg")
    set(FLOW_BUILD_TYPE debug)
    MESSAGE(STATUS "Auto set FLOW_BUILD_TYPE to ${FLOW_BUILD_TYPE}")
  elseif ("${DOCKER_VERSION}" MATCHES "rel")
    set(FLOW_BUILD_TYPE release)
    MESSAGE(STATUS "Auto set FLOW_BUILD_TYPE to ${FLOW_BUILD_TYPE}")
  else()
    MESSAGE(FATAL_ERROR "Could not decide on FLOW_BUILD_TYPE from dockerversion '${DOCKER_VERSION}', need configuration file: ${CONFIG_CMAKE}")
  endif()
endif()

# by default reuse number of cpus for make from environment
if (NOT MAKE_NUMCPUS)
    set(MAKE_NUMCPUS $ENV{MAKE_NUMCPUS})
endif()

if (MAKE_NUMCPUS)
    set(MAKE_NUMCPUS_FLAG "-j ${MAKE_NUMCPUS}")
endif()    



# If the CCACHE_CXX_COMPILER is set 
# then the link to this compiler with ccache is created and used.
if(CCACHE_CXX_COMPILER)
    set (CMAKE_CXX_COMPILER "ccache_c++")
    configure_file(${CMAKE_SOURCE_DIR}/CMake/ccache_c++_link_template ${CMAKE_BINARY_DIR}/ccache_c++)
endif(CCACHE_CXX_COMPILER)


# try to use ccache for compiling and linking if flag is set
if(USE_CCACHE)
    find_program(CCACHE_FOUND ccache)
    if(CCACHE_FOUND)
        set_property(GLOBAL PROPERTY RULE_LAUNCH_COMPILE ccache)
        set_property(GLOBAL PROPERTY RULE_LAUNCH_LINK ccache)
    endif(CCACHE_FOUND)
endif(USE_CCACHE)


# Project name and language(s)
project(Flow123d CXX C Fortran)

# CMake version required for pre-build this project
cmake_minimum_required(VERSION 2.8.8)
    
cmake_policy(SET CMP0017 NEW)
if (${CMAKE_VERSION} VERSION_GREATER 3.0)
    cmake_policy(SET CMP0037 OLD)
endif()


# allow duplicate custom targets for makefiles
set_property(GLOBAL PROPERTY ALLOW_DUPLICATE_CUSTOM_TARGETS 1)
set(ALLOW_DUPLICATE_CUSTOM_TARGETS 1) #variable that can be tested


################################################################################# set CMake Modules
# Set all module paths
set(CMAKE_MODULE_PATH "${Flow123d_SOURCE_DIR}/CMake/Modules" "${CMAKE_ROOT}/Modules")

if (NOT EXTERNAL_PROJECT_DIR)
    # default directory for building external libraries 
    set(EXTERNAL_PROJECT_DIR ${PROJECT_SOURCE_DIR}/../external_projects)
endif()


############################################################################################# set compiler
# check C++11 properties of the compiler
# Initialize CXXFLAGS.




##############
# set appropriate compiler flags for debug/release compilation modes
message(STATUS "Flow build type: ${FLOW_BUILD_TYPE} C++ compiler falgs: ${FLOW_CC_FLAGS}")
include (Flow123d_config_header)




  
message(STATUS "Using build type: ${FLOW_BUILD_TYPE}")
if(FLOW_BUILD_TYPE STREQUAL "debug")
  # debug mode
  set(CMAKE_CXX_FLAGS "-g -O0")
  set(CMAKE_C_FLAGS "-g -O0")
  set(CMAKE_EXE_LINKER_FLAGS "-g ")
  set(CMAKE_BUILD_TYPE Debug)
  flow_define(DEBUG)
elseif(FLOW_BUILD_TYPE STREQUAL "release")
  # Release mode
  set(CMAKE_CXX_FLAGS "-O3 -DNDEBUG" )
  set(CMAKE_C_FLAGS "-O3 -DNDEBUG")
  set(CMAKE_BUILD_TYPE Release)
  # armadillo is configured using this definition
  add_definitions(-DARMA_NO_DEBUG)    # no bound checks in armadillo 
  flow_define(DEBUG_PROFILER)
  flow_define(RUN_UNIT_BENCHMARKS)    # build benchmark unit tests
elseif(FLOW_BUILD_TYPE STREQUAL "release_with_debug")
  # Release mode with debug info 
  set(CMAKE_CXX_FLAGS "-O3")
  set(CMAKE_C_FLAGS "-O3")
  set(CMAKE_BUILD_TYPE Debug)
  flow_define(DEBUG)
elseif(FLOW_BUILD_TYPE STREQUAL "coverage")
  # Coverage mode creates instrumented binary
  SET(CMAKE_CXX_FLAGS "-g -O0 --coverage")
  SET(CMAKE_C_FLAGS "-g -O0 --coverage")
  SET(CMAKE_EXE_LINKER_FLAGS "-g -O0 --coverage")
  set(CMAKE_BUILD_TYPE Debug)
  
  # set external libs flags
  SET(EXTERNAL_LIBS_CXX_FLAGS "-g -O0")
  SET(EXTERNAL_LIBS_C_FLAGS "-g -O0")
  
  flow_define(DEBUG)
  # INCLUDE(CodeCoverage)
  # add custom target to collect coverage into xml file
  # COVERAGE_SETUP_COBERTURA_COLLECTION (cobertura-collect-coverage coverage.xml)
else()
endif()

if(FLOW_CC_FLAGS)
  # user given flags overrides everything  
  message(STATUS "Using user specified compiler options: ${FLOW_CC_FLAGS}")
  set(CMAKE_CXX_FLAGS  "${FLOW_CC_FLAGS}")
  set(CMAKE_C_FLAGS  "${FLOW_CC_FLAGS}")
  set(CMAKE_EXE_LINKER_FLAGS "${FLOW_CC_FLAGS}")
endif(FLOW_CC_FLAGS)

# unset predefined build type resolution variables to have better control over
# used compiler flags
set(CMAKE_C_FLAGS_${CMAKE_BUILD_TYPE} "")
set(CMAKE_CXX_FLAGS_${CMAKE_BUILD_TYPE} "")
set(CMAKE_FORTRAN_FLAGS_${CMAKE_BUILD_TYPE} "")

# Compiler-specific C++14 activation.
set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++14")
      # -Wvla - variable length arrays
      set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -pedantic -Wall -Wextra -Wno-vla -Wno-unused-local-typedefs -Wno-deprecated-declarations")

## Compiler-specific C++11 activation.
## It was replaced in commit a66193035 with C++14 specification.
## Needs special definition of trigonometric functions held in the comment below:
# if ("${CMAKE_CXX_COMPILER_ID}" MATCHES "GNU")        
#     ...
#     # Necessary for FParser - TODO: move to better place
#     # ?? what about Clang
#     flow_define(HAVE_SINCOS)
# else()
#     ...
# endif()

# set encoding of VTU output files to Little Endian, VTU files can be encoded in Little Endian or Big Endian format  
flow_define(LITTLE_ENDIAN true)


#########################################################################
# Determining timer with highest resolution
include(ProfilerTimerCheck)
set(TIMER_RESOLUTION No)

# try to get best timer on windows
if (WIN32 OR CYGWIN)
    # on windows and cygwin try to get access QueryPerformanceCounter timer
    TEST_TIMER (HAVE_TIMER_QUERY_PERFORMANCE_COUNTER
      ${CMAKE_SOURCE_DIR}/CMake/TimerTest/timer_queryperformance.cpp) 

    # if success save theese values via flow_define and flow_define_constant macros
    if (HAVE_TIMER_QUERY_PERFORMANCE_COUNTER AND TIMER_RESOLUTION)
        flow_define           (HAVE_TIMER_QUERY_PERFORMANCE_COUNTER)
        # flow_define_constant  (TIMER_RESOLUTION "${TIMER_RESOLUTION}")

        message(STATUS "Timer used:       QueryPerformanceCounter")
        message(STATUS "Timer resolution: ${TIMER_RESOLUTION} [us]")
    else()
        MESSAGE (WARNING "Couldn't access QueryPerformanceCounter timer")   
    endif()
endif()

if (NOT TIMER_RESOLUTION)
    # otherwise try to access chrono high resolution timer (multiplatform)
    # expecting c++11 support
    TEST_TIMER (HAVE_TIMER_CHRONO_HIGH_RESOLUTION
      ${CMAKE_SOURCE_DIR}/CMake/TimerTest/timer_chrono_highresolution.cpp)

    # if success save theese values via flow_define and flow_define_constant macros
    if (HAVE_TIMER_CHRONO_HIGH_RESOLUTION AND TIMER_RESOLUTION)
        flow_define           (HAVE_TIMER_CHRONO_HIGH_RESOLUTION)
        # flow_define_constant  (TIMER_RESOLUTION "${TIMER_RESOLUTION}")

        message(STATUS "Timer used:       chrono highresolution timer")
        message(STATUS "Timer resolution: ${TIMER_RESOLUTION} [us]")
    else()
        MESSAGE (FATAL_ERROR "Couldn't access chrono timer")    
    endif()
endif ()



#########################################################################
# Check that glibc has working execinfo.h
include(CheckIncludeFile)
check_include_file("execinfo.h" EXEC_INFO_H)
if (EXEC_INFO_H) 
    flow_define(HAVE_EXEC_INFO)
    
    MESSAGE(STATUS "Looking for abi::__cxa_demangle")
    TRY_COMPILE(HAVE_DEMAGLER
        ${CMAKE_BINARY_DIR}
        ${CMAKE_SOURCE_DIR}/CMake/CheckDemagler.cpp
    )
    
    IF(HAVE_DEMAGLER)
      MESSAGE(STATUS "Looking for abi::__cxa_demangle - found")
      flow_define(HAVE_DEMAGLER)
    ELSE(HAVE_DEMAGLER)
      MESSAGE(STATUS "Looking for abi::__cxa_demangle - not found")
    ENDIF(HAVE_DEMAGLER)
    
endif(EXEC_INFO_H)



# set value for external libs flags if specified, otherwise use default cmake flags
INCLUDE (MacroValidValue)
SET_VALID_VALUE ("EXTERNAL_LIBS_C_FLAGS"         ${EXTERNAL_LIBS_C_FLAGS}         ${CMAKE_C_FLAGS})
SET_VALID_VALUE ("EXTERNAL_LIBS_CXX_FLAGS"       ${EXTERNAL_LIBS_CXX_FLAGS}       ${CMAKE_CXX_FLAGS})
SET_VALID_VALUE ("EXTERNAL_LIBS_Fortran_FLAGS"   ${EXTERNAL_LIBS_Fortran_FLAGS}   ${CMAKE_Fortran_FLAGS})


# store compiler flags (into config.h file)
flow_define(COMPILER_FLAGS_ ${CMAKE_CXX_FLAGS})
# store current source dir (for profiler purposes)
flow_define(SOURCE_DIR ${CMAKE_SOURCE_DIR})

set(JSON_BuildTests OFF CACHE INTERNAL "")
add_subdirectory("${CMAKE_SOURCE_DIR}/third_party/json-3.7.1")

######################################################################### find libraries
#
# find boost or use our own
#
# - Try to find Boost include dirs and libraries
#
# NOTE: Take note of the Boost_ADDITIONAL_VERSIONS variable below.
# Due to Boost naming conventions and limitations in CMake this find
# module is NOT future safe with respect to Boost version numbers,
# and may break.
#
# Input variables:
#
#   Boost_USE_MULTITHREADED      Can be set to OFF to use the non-multithreaded
#   Boost_USE_STATIC_LIBS        Can be set to ON to force the use of the static
#   Boost_NO_SYSTEM_PATHS        Set to TRUE to suppress searching in system (or other locations outside of BOOST_ROOT or BOOST_INCLUDEDIR).
#   Boost_USE_STATIC_RUNTIME     If enabled, searches for boost libraries linked against a static C++ standard library
#   Boost_USE_DEBUG_PYTHON       
#   Boost_USE_STLPORT            If enabled, searches for boost libraries compiled against the STLPort standard
#   Boost_USE_STLPORT_DEPRECATED_NATIVE_IOSTREAMS    If enabled, searches for boost libraries compiled against the deprecated STLPort "native iostreams" feature ('n' ABI tag).
#   Boost_ADDITIONAL_VERSIONS    A list of version numbers to use for searching
#   Boost_DEBUG                  Set this to TRUE to enable debugging output
#   Boost_DETAILED_FAILURE_MSG   FindBoost doesn't output detailed information
#   Boost_COMPILER               Set this to the compiler suffix used by Boost
#   Boost_THREADAPI                When building boost.thread, sometimes the name of the library contains an additional "pthread" or "win32"
#   Boost_REALPATH               Resolves symbolic links for discovered boost libraries
#
#
#   BOOST_ROOT or BOOSTROOT      The preferred installation prefix for searching for
#                                back on the system paths, set Boost_NO_SYSTEM_PATHS
#   BOOST_INCLUDEDIR             Set this to the include directory of Boost, if the
#   BOOST_LIBRARYDIR             Set this to the lib directory of Boost, if the
#
# Variables defined by this module:
#
#   Boost_FOUND                         System has Boost, this means the include dir was
#   Boost_INCLUDE_DIRS                  Boost include directories: not cached
#   Boost_INCLUDE_DIR                   This is almost the same as above, but this one is cached and may be modified by advanced users
#   Boost_LIBRARIES                     Link to these to use the Boost libraries that you
#   Boost_LIBRARY_DIRS                  The path to where the Boost library files are.
#   Boost_VERSION                       The version number of the boost libraries that
#   Boost_LIB_VERSION                   The version number in filename form as
#   Boost_MAJOR_VERSION                 major version number of boost
#   Boost_MINOR_VERSION                 minor version number of boost
#   Boost_SUBMINOR_VERSION              subminor version number of boost
#   Boost_LIB_DIAGNOSTIC_DEFINITIONS    [WIN32 Only] You can call
#   Boost_${COMPONENT}_FOUND            True IF the Boost library "component" was found.
#   Boost_${COMPONENT}_LIBRARY          Contains the libraries for the specified Boost

message(STATUS "=======================================================")
message(STATUS "====== BOOST ==========================================")
message(STATUS "=======================================================")
set(Boost_USE_STATIC_LIBS       ON)
set(Boost_USE_MULTITHREADED     OFF)
set(Boost_USE_STATIC_RUNTIME    OFF)
set(EXTERNAL_BOOST_DIR  "${EXTERNAL_PROJECT_DIR}/boost_build")  
#set(REQUIRED_Boost_COMPONENTS program_options serialization regex filesystem system)
set(REQUIRED_Boost_COMPONENTS program_options filesystem) 
 
# set(Boost_DEBUG TRUE)
# For some mysterious reason, the second find_package fails to find
# components of boost in our own build, if we hint BOOST_ROOT of our build to the
# first find_package.
#
# The root of the error is in FindBoost, when it calls find_library for individual components.
# Success scenario:
# 1) First find_package didn't found Boost (in system direcotries)
# 2) we install our Boost
# 3) the second find_package calls find_library to find boost components and succed
#
# Fail scenario:
# 1) First find_package didn't found Boost (we set BOOST_ROOT to the empty target directory of the external project)
# 2) we install our Boost
# 3) the second find_package calls find_library to find boost components and fail
#    the parameters of the find_library are exactly same as in the previous case.
#

# If no explicit setting from config.cmake ....
if (NOT BOOST_ROOT)
    set(BOOST_ROOT $ENV{BOOST_ROOT})
endif()    

# If no system setting ...
if (NOT BOOST_ROOT)
	set(BOOST_ROOT ${EXTERNAL_BOOST_DIR})
endif()

if (NOT Boost_FORCE_REBUILD)
    find_package(Boost 1.40.0 COMPONENTS ${REQUIRED_Boost_COMPONENTS} )
endif()    

if (NOT Boost_FOUND) 
  message(STATUS "Try to install our own Boost.")
  message("NCPU: ${MAKE_NUMCPUS}")
  include(InstallBoost)  
  
  # we try to erease from the cache everything that can cause failure of the second find_package
  # however just unsetting Boost_INCLUDE_DIR should be enough
  
 
  unset(Boost_FOUND CACHE) 
  unset(Boost_INCLUDE_DIR CACHE)        # overwrite cache    
  unset(Boost_LIBRARY_DIRS CACHE)
  unset(Boost_LIB_VERSION CACHE)
  unset(Boost_VERSION CACHE)
  foreach(COMPONENT ${REQUIRED_Boost_COMPONENTS})
      string(TOUPPER ${COMPONENT} COMPONENT)

      unset(Boost_${COMPONENT}_FOUND CACHE)
      unset(Boost_${COMPONENT}_LIBRARY CACHE)
      unset(Boost_${COMPONENT}_LIBRARY_DEBUG CACHE)
      unset(Boost_${COMPONENT}_LIBRARY_RELEASE CACHE)
            
  endforeach(COMPONENT)
  
  set(BOOST_ROOT ${EXTERNAL_BOOST_DIR})    
  find_package( Boost 1.42.0 REQUIRED COMPONENTS ${REQUIRED_Boost_COMPONENTS}  )
   
endif()

message(STATUS "-------------------------------------------------------")
message(STATUS "Boost_FOUND = ${Boost_FOUND}")
message(STATUS "BOOST_ROOT = ${BOOST_ROOT}")
message(STATUS "Boost_LIBRARIES = ${Boost_LIBRARIES}")
message(STATUS "Boost_LIBRARY_DIRS = ${Boost_LIBRARY_DIRS}")
message(STATUS "Boost_INCLUDE_DIR = ${Boost_INCLUDE_DIR}")
message(STATUS "=======================================================\n\n")


# We are looking for PETSC, f succed this module set:
#  PETSC_FOUND        - system has PETSc
#  PETSC_INCLUDES     - the PETSc include directories
#  PETSC_LIBRARIES    - Link these to use PETSc
#  PETSC_COMPILER     - Compiler used by PETSc, helpful to find a compatible MPI
#  PETSC_DEFINITIONS  - Compiler switches for using PETSc
#  PETSC_MPIEXEC      - Executable for running MPI programs
#  PETSC_VERSION      - Version string (MAJOR.MINOR.SUBMINOR)
message(STATUS "=======================================================")
message(STATUS "====== PETSC ==========================================")
message(STATUS "=======================================================")
include (FindPackageMultipass)

# this export list is for possible BDDCML installation
set(PETSC_EXPORT_LIST 
  MPIFC_SHOW MPICXX_SHOW MPICC_SHOW CC FC CC_FLAGS PCC_FLAGS FC_FLAGS CPP_FLAGS
  METIS_INCLUDE BLASLAPACK_INCLUDE PARMETIS_INCLUDE MUMPS_INCLUDE
  MPI_LIB MUMPS_LIB  METIS_LIB PARMETIS_LIB BLASLAPACK_LIB SCALAPACK_LIB
  UMFPACK_LIB)  

message(STATUS "PETSC_DIR = ${PETSC_DIR}")
message(STATUS "PETSC_ARCH = ${PETSC_ARCH}")

find_package(PETSc 3.6.0 REQUIRED)
flow_define(HAVE_PETSC)
flow_define(HAVE_MPI)

set(MPI_CXX_LIBRARIES ${PETSC_EXTERNAL_LIB})
set(MPI_CXX_INCLUDE_PATH ${PETSC_INCLUDES})

message(STATUS "-------------------------------------------------------")
message(STATUS "PETSC_FOUND = ${PETSC_FOUND}")
message(STATUS "PETSC_LIBS = ${PETSC_LIBRARIES}")
message(STATUS "PETSC_EXTERNAL_LIB = ${PETSC_EXTERNAL_LIB}")
message(STATUS "=======================================================\n\n")


################################################################################
#  BDDCML_FOUND - set to true if the library is found
message(STATUS "=======================================================")
message(STATUS "====== BDDCML =========================================")
message(STATUS "=======================================================")
message(STATUS "BDDCML_ROOT = ${BDDCML_ROOT}")

find_package(BDDCML)
if (BDDCML_FOUND)
    flow_define(HAVE_BDDCML)
    message(STATUS "-------------------------------------------------------")
    message(STATUS "BDDCML_FOUND = ${BDDCML_FOUND}")
    message(STATUS "BDDCML_INCLUDES = ${BDDCML_INCLUDES}")
    message(STATUS "BDDCML_LIBRARY = ${BDDCML_LIBRARY}")
    message(STATUS "=======================================================\n\n") 
endif()



#################################################################################
#  Armadillo_FOUND - set to true if the library is found
#  Armadillo_INCLUDE_DIRS - list of required include directories
#  Armadillo_LIBRARIES - list of libraries to be linked
#  Armadillo_CONFIGURE_CMAKE - path to configure cmake script that introduce target 'armadillo'
#                              with correctly set dependencies on other libraries
message(STATUS "=======================================================")
message(STATUS "====== ARMADILLO ======================================")
message(STATUS "=======================================================")
message(STATUS "Armadillo_ROOT_HINT = ${Armadillo_ROOT_HINT}")


find_package(ArmadilloFlow 3.2.4 REQUIRED)
include("${Armadillo_CONFIGURE_FILE}")   # this introduce target "armadillo" with correct dependencies
    # until we know how to take only blas and lapack from PETSC we have to make Armadillo depend on whole PETSC lib.
    # this may be fragile because of _NOCONFIG
set_target_properties(armadillo PROPERTIES
  IMPORTED_LINK_INTERFACE_LIBRARIES_NOCONFIG "${PETSC_EXTERNAL_LIB}")

message(STATUS "-------------------------------------------------------")
message(STATUS "Armadillo_FOUND = ${Armadillo_FOUND}")
message(STATUS "Armadillo_LIBRARIES = ${Armadillo_LIBRARIES}")
message(STATUS "Armadillo_INCLUDE_DIRS = ${Armadillo_INCLUDE_DIRS}")
message(STATUS "Armadillo_CONFIGURE_FILE = ${Armadillo_CONFIGURE_FILE}")
message(STATUS "=======================================================\n\n")


#################################################################################
#  YAML_CPP_FOUND - set to true if the library is found
message(STATUS "=======================================================")
message(STATUS "====== YAML-CPP =======================================")
message(STATUS "=======================================================")
message(STATUS "YamlCpp_ROOT_HINT = ${YamlCpp_ROOT_HINT}")

find_package(YamlCpp REQUIRED)

message(STATUS "-------------------------------------------------------")
message(STATUS "YamlCpp_FOUND = ${YamlCpp_FOUND}")
message(STATUS "YamlCpp_LIBRARY = ${YamlCpp_LIBRARY}")
message(STATUS "YamlCpp_INCLUDE_DIR = ${YamlCpp_INCLUDE_DIR}")
message(STATUS "=======================================================\n\n")


#################################################################################
#  ZLIB_FOUND - set to true if the library is found
message(STATUS "=======================================================")
message(STATUS "====== ZLIB-CORE ======================================")
message(STATUS "=======================================================")
message(STATUS "Zlib_ROOT_HINT = ${Zlib_ROOT_HINT}")

find_package(Zlib)

# use zlib compression of VTK files
if(Zlib_FOUND)
    flow_define(HAVE_ZLIB)
endif()

message(STATUS "-------------------------------------------------------")
message(STATUS "Zlib_FOUND = ${Zlib_FOUND}")
message(STATUS "Zlib_LIBRARY = ${Zlib_LIBRARY}")
message(STATUS "Zlib_INCLUDE_DIR = ${Zlib_INCLUDE_DIR}")
message(STATUS "=======================================================\n\n")


####################################################################################
# PYTHON
message(STATUS "=======================================================")
message(STATUS "====== PYTHON =========================================")
message(STATUS "=======================================================")

# set PYTHON_EXTRA_MODULES_PATH to lib/python/flow123d
# if we are in installed flow123d, this direcoty exists and 
# contains all python flow123d files
# otherwise
# we create link to src so python_loader is able to use modules
set(PYTHON_EXTRA_MODULES_PATH lib/python/flow123d)
# PYTHONPATH is relative to cwd
# need to make changes in c++ to correctly alter PYTHONPATH relative to binary location
# for now, we also add absolute location
# now it contains:
#       relative path
#       absolute path
#       /opt/flow123d - for production
#       /opt/flow123d/flow123d - for packaging
flow_define(PYTHON_EXTRA_MODULES_PATH ../${PYTHON_EXTRA_MODULES_PATH}:${CMAKE_SOURCE_DIR}/${PYTHON_EXTRA_MODULES_PATH}:/opt/flow123d/${PYTHON_EXTRA_MODULES_PATH}:/opt/flow123d/flow123d/${PYTHON_EXTRA_MODULES_PATH})
execute_process(COMMAND mkdir -p ${CMAKE_SOURCE_DIR}/lib/python)
if (COPY_INSTEAD_OF_SYMLINK)
	execute_process(COMMAND
	    ${CMAKE_COMMAND} -E copy ${CMAKE_SOURCE_DIR}/src/python ${CMAKE_SOURCE_DIR}/lib/python/flow123d
	    WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
	)
	execute_process(COMMAND
	    ${CMAKE_COMMAND} -E copy_directory ${CMAKE_SOURCE_DIR}/src/python ${CMAKE_SOURCE_DIR}/lib/python/flow123d
	    WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
	)
	
	set(EXECUTABLE_PATH ${MPIEXEC_PATH})
configure_file(${CMAKE_SOURCE_DIR}/CMake/executable_link_template ${CMAKE_BINARY_DIR}/bin/mpiexec)
	
else()
	execute_process(COMMAND
	    ${CMAKE_COMMAND} -E create_symlink ${CMAKE_SOURCE_DIR}/src/python ${CMAKE_SOURCE_DIR}/lib/python/flow123d
	    WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
	)
endif()

# use python on demand
if(USE_PYTHON)
    flow_define(HAVE_PYTHON)
endif()

# find python 
set(Python_ADDITIONAL_VERSIONS 3.4) # for cmake version 2.8
find_package(PythonInterp 3.4 REQUIRED)
find_package(PythonLibs   3.4 REQUIRED)
# TODO: remove dependency on python
include_directories(${PYTHON_INCLUDE_DIRS})

# create shortcuts for runtest

# use relative path
set(PYTHON_RUNTEST ../bin/python/runtest.py)
#configure_file(${CMAKE_SOURCE_DIR}/CMake/unix_runtest_template ${CMAKE_SOURCE_DIR}/bin/runtest @ONLY)
#configure_file(${CMAKE_SOURCE_DIR}/CMake/unix_runtest_template ${CMAKE_SOURCE_DIR}/tests/runtest @ONLY)
#configure_file(${CMAKE_SOURCE_DIR}/CMake/unix_runtest_template ${CMAKE_SOURCE_DIR}/tests_old/runtest @ONLY)

# make shortcut also for exec_parallel.py
set(PYTHON_RUNTEST ../bin/python/exec_parallel.py)
configure_file(${CMAKE_SOURCE_DIR}/CMake/unix_runtest_template ${CMAKE_SOURCE_DIR}/bin/exec_parallel @ONLY)

# restore absolute path
set(PYTHON_RUNTEST ${CMAKE_SOURCE_DIR}/bin/python/runtest.py)

message(STATUS "-------------------------------------------------------")
message(STATUS "PYTHONINTERP_FOUND = ${PYTHONINTERP_FOUND}")
message(STATUS "PYTHON_EXECUTABLE = ${PYTHON_EXECUTABLE}")
message(STATUS "PYTHON_VERSION_STRING = ${PYTHON_VERSION_STRING}")
message(STATUS "-------------------------------------------------------")
message(STATUS "PYTHONLIBS_FOUND = ${PYTHONLIBS_FOUND}")
message(STATUS "PYTHON_LIBRARIES = ${PYTHON_LIBRARIES}")
message(STATUS "PYTHON_INCLUDE_DIRS = ${PYTHON_INCLUDE_DIRS}")
message(STATUS "PYTHONLIBS_VERSION_STRING = ${PYTHONLIBS_VERSION_STRING}")
message(STATUS "=======================================================\n\n")


#################################################################################
#  PUGI_XML_FOUND - set to true if the library is found
message(STATUS "=======================================================")
message(STATUS "====== PUGI-XML =======================================")
message(STATUS "=======================================================")
message(STATUS "PugiXml_ROOT_HINT = ${PugiXml_ROOT_HINT}")

find_package(PugiXml REQUIRED)

message(STATUS "-------------------------------------------------------")
message(STATUS "PugiXml_FOUND = ${PugiXml_FOUND}")
message(STATUS "PugiXml_LIBRARY = ${PugiXml_LIBRARY}")
message(STATUS "PugiXml_INCLUDE_DIR = ${PugiXml_INCLUDE_DIR}")
message(STATUS "=======================================================\n\n")


message(STATUS "Installing HTML documentation")
include(InstallHTMLDOC)


message(STATUS "\n ======================================================= ")


# directory containing includes  (header files)
set(CONFIG_FILE_INCLUDES_DIRECTORY "${CMAKE_BINARY_DIR}/includes")

######################################################################################################################### 
# set includes and libraries paths

include_directories(
    ${CMAKE_CURRENT_BINARY_DIR}                            # include the output directory, where the rev_num.h file is generated    
    ${CMAKE_SOURCE_DIR}/src
    ${CONFIG_FILE_INCLUDES_DIRECTORY}
    #${CMAKE_SOURCE_DIR}/third_party/metis-include
    ${CMAKE_SOURCE_DIR}/third_party/fparser-4.5.1
    ${CMAKE_SOURCE_DIR}/third_party/FADBAD++
    ${PETSC_INCLUDES}
    ${PYTHON_INCLUDE_DIRS}
    
    ${Boost_INCLUDE_DIR}
    ${Armadillo_INCLUDE_DIRS}
    ${MPI_CXX_INCLUDE_PATH}
    ${YamlCpp_INCLUDE_DIR}
    ${PugiXml_INCLUDE_DIR}
    ${Zlib_INCLUDE_DIR}
    ${CMAKE_BINARY_DIR}/src/dealii/include    # deal generates config.h
    ${CMAKE_SOURCE_DIR}/src/dealii/include
#    ${CMAKE_SOURCE_DIR}/third_party/tbb43_20150316oss/include
)

if (BDDCML_FOUND)
    include_directories(${BDDCML_INCLUDES})
endif()    
    
######################################################################################################################
# custom targets


#### Safety check of build_tree link
# This check is called by every target, if it fails we are probably running makefile generated for another branch.
# So we must break the build.
add_custom_target(update_source_tree 
    COMMAND ${CMAKE_SOURCE_DIR}/bin/git_post_checkout_hook
    WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
)

#### create rev_num.h - file with some info about code origin
# a custom target - is always built -> call cmake script to generate rev_num.h
# Get the current working branch

# This command is run from build dir and cmake set CMAKE_SOURCE_DIR to the build dir; thus we have to 
# pass different variable with source dir, namely FLOW123D_SOURCE_DIR
set(CONFIG_FILE_REV_NUM "${CONFIG_FILE_INCLUDES_DIRECTORY}/rev_num.h")
add_custom_target(rev_num_always 
    COMMAND ${CMAKE_COMMAND}  -DFLOW123D_SOURCE_DIR=${CMAKE_SOURCE_DIR} -DOUTPUT_FILE_PATH=${CONFIG_FILE_REV_NUM} -P ${CMAKE_SOURCE_DIR}/CMake/Modules/Flow123d_git_rev_num.cmake 
)
# rev_num.h is a generated file
set_source_files_properties(${CONFIG_FILE_REV_NUM}
    PROPERTIES GENERATED TRUE
    HEADER_FILE_ONLY TRUE)

### create config.h - file containing definitions flags and other information
# called everytime
set(CONFIG_FILE_HEADER "${CONFIG_FILE_INCLUDES_DIRECTORY}/config.h")
GENERATE_CONFIG_H ("${CONFIG_FILE_HEADER}")
# config.h is a generated file
set_source_files_properties(${CONFIG_FILE_HEADER}
    PROPERTIES GENERATED TRUE
    HEADER_FILE_ONLY TRUE)


# get git information also through configure phase
set(FLOW123D_SOURCE_DIR ${CMAKE_SOURCE_DIR})
include(Flow123d_git_info)

# version in reference manual
STRING(REGEX REPLACE "_" "\\\\_" GIT_VERSION_FULL_LATEX ${GIT_VERSION_FULL})
configure_file(${CMAKE_SOURCE_DIR}/CMake/flow_version_tex_template ${CMAKE_SOURCE_DIR}/doc/reference_manual/flow_version.tex)    
    
#### make a script which call correct mpiexec (taken from PETSC)
# we can not use link, since mpiexec can be sript which reads its location from calling sequence
#

#message(STATUS "petsc_mpiexec: ${PETSC_MPIEXEC}")
if(PETSC_MPIEXEC)     
      # missing CMAKE command to test working executable, moreover it could be without full path
      
      # resolve symlinks
      get_filename_component(RESOLVED_MPIEXEC "${PETSC_MPIEXEC}" REALPATH)
      
      # we verify that REALPATH of PETSC_MPIEXEC exists
      # if not, PETSC_MPIEXEC is probably only command 'mpiexec'
      # in that case we use this value, since PETSc was ok with that
      if(EXISTS RESOLVED_MPIEXEC)
          message(STATUS "using petsc_mpiexec: ${RESOLVED_MPIEXEC}")
          set(MPIEXEC_PATH "${RESOLVED_MPIEXEC}")
      else()
          message(STATUS "path does not point to valid location: ${RESOLVED_MPIEXEC} using ${PETSC_MPIEXEC}")
          # set MPIEXEC_PATH to full path of PETSC_MPIEXEC
          execute_process(COMMAND
              which ${PETSC_MPIEXEC}
              OUTPUT_VARIABLE MPIEXEC_PATH OUTPUT_STRIP_TRAILING_WHITESPACE
          )
          get_filename_component(MPIEXEC_PATH "${MPIEXEC_PATH}" REALPATH)
      endif()
      
else()
      find_program(SYSTEM_MPIEXEC mpiexec)
      if(SYSTEM_MPIEXEC)
        message(WARNING "Missing mpiexec in PETSc instalation. Using system wide mpiexec.") 
        get_filename_component(RESOLVED_MPIEXEC "${SYSTEM_MPIEXEC}" REALPATH)
        set(MPIEXEC_PATH "${RESOLVED_MPIEXEC}")
      else()
        message(WARNING "Missing any mpiexec.")
      endif()
endif()
# verify that mpiexec exists (on CYGWIN we have to consider .exe suffix)
if( CYGWIN AND NOT EXISTS ${MPIEXEC_PATH})
    if( EXISTS "${MPIEXEC_PATH}.exe")
        set(MPIEXEC_PATH "${MPIEXEC_PATH}.exe")
    endif()
endif()

# create shell script which call correct mpiexec binary
set(EXECUTABLE_PATH ${MPIEXEC_PATH})
configure_file(${CMAKE_SOURCE_DIR}/CMake/executable_link_template ${CMAKE_BINARY_DIR}/bin/mpiexec)

##################################################################################################################### 
# Macro to make usefull links in SOURCE_DIR
##################################


# initialize a global variable (cached)
set(SOURCE_TREE_LINKS  "" CACHE INTERNAL "" FORCE)
              

macro(make_link_in_source FILE)
  # get relative filename if absolute address is given
  string(REPLACE "${CMAKE_SOURCE_DIR}/" "" rel_path "${FILE}")
  if (COPY_INSTEAD_OF_SYMLINK)
    execute_process(
      COMMAND ${CMAKE_COMMAND} -E copy ${CMAKE_SOURCE_DIR}/build_tree/${rel_path} ${CMAKE_SOURCE_DIR}/${rel_path}
      WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
    )
    execute_process(
      COMMAND ${CMAKE_COMMAND} -E copy_directory ${CMAKE_SOURCE_DIR}/build_tree/${rel_path} ${CMAKE_SOURCE_DIR}/${rel_path}
      WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
    )
  else()
    execute_process(
      COMMAND ${CMAKE_COMMAND} -E create_symlink ${CMAKE_SOURCE_DIR}/build_tree/${rel_path} ${CMAKE_SOURCE_DIR}/${rel_path}
      WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
    )
  endif()
  
  set(SOURCE_TREE_LINKS  ${SOURCE_TREE_LINKS} ${CMAKE_SOURCE_DIR}/${FILE} CACHE INTERNAL "")
endmacro()

# use wreappers instead of links (or copy) for executables
#set(EXECUTABLE_PATH ${CMAKE_BINARY_DIR}/bin/flow123d)
#configure_file(${CMAKE_SOURCE_DIR}/CMake/executable_link_template ${CMAKE_SOURCE_DIR}/bin/flow123d)
#set(EXECUTABLE_PATH ${CMAKE_BINARY_DIR}/bin/mpiexec)
#configure_file(${CMAKE_SOURCE_DIR}/CMake/executable_link_template ${CMAKE_SOURCE_DIR}/bin/mpiexec)

#make_link_in_source(bin/flow123d)
#make_link_in_source(bin/mpiexec)

make_link_in_source(bin/tests/Makefile)
make_link_in_source(unit_tests/Makefile)

############################################################33 process subdirectories


# add fparser sources
message(STATUS "add subdir FParser")
add_subdirectory(third_party/fparser-4.5.1)

# add tbb sources
#message(STATUS "add subdir TBB")
#add_subdirectory(third_party/tbb43_20150316oss)

# add various our sources
add_subdirectory(src)
add_subdirectory(doc)
add_subdirectory(unit_tests) # this also adds third_party/gtest

# install benchmark miniapp
execute_process(COMMAND
    make -C ${CMAKE_SOURCE_DIR}/src/python/benchmarks all
)

# optionally make targets for building reference_manual
# we have to call separate instance of CMake since UseLatex.cmake
# depends on many tools. If any of them miss CMake fails, while we can live
# without targets for Latex doc.
message(STATUS "==== child CMake for Latex sources ====") 
execute_process(COMMAND ${CMAKE_COMMAND} -E make_directory ${CMAKE_BINARY_DIR}/doc/reference_manual)
#execute_process(COMMAND ${CMAKE_COMMAND} ${CMAKE_SOURCE_DIR}/doc/reference_manual
#                WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/doc/reference_manual
#                RESULT_VARIABLE res_var)
message(STATUS "==== back in main CMake ====")
#if(NOT "${res_var}" STREQUAL "0")
#  message(WARNING "Some problem with Latex, 'make ref-doc' may not work.")
#endif()


#################################################################### add targets


if( PETSC_FOUND )
  ########################
# Setting for RPATH

 
  SET (CMAKE_EXE_LINKER_FLAGS
    "${CMAKE_EXE_LINKER_FLAGS} -Wl,-rpath,$ORIGIN/../lib" )
  
  # bug in cmake https://cmake.org/Bug/view.php?id=13934
  # causing error during CPack operation
  # 
  #    CMake Error at cmake_install.cmake (FILE):
  #        file INSTALL cannot find
  #        
  # This is workaround posted on cmake.org. This bug occured in version 2.8 
  # (first occurence in version 2.8.10) but it is platform depended, so some 
  # unix systems will not have this bug
  # 
  # By setting CMAKE_EXECUTABLE_FORMAT to "ELF" in the toolchain file you 
  #   tell CMake that it can use its builtin ELF editor instead of relinking
  if (${CMAKE_VERSION} VERSION_GREATER 2.8.10 AND ${CMAKE_VERSION} VERSION_LESS 2.9)
      set(CMAKE_EXECUTABLE_FORMAT "ELF")
  endif()
  
  ### add target for flow123d binary
  add_executable(bin/flow123d src/main.cc)
  target_link_libraries(bin/flow123d 
        flow123d_lib 
        ${PETSC_ADDITIONAL_LIBS}
        ${PYTHON_LIBRARIES}
        ${Boost_LIBRARIES}
  )     # set libraries linked to main binary (order matter)
  add_dependencies(bin/flow123d update_source_tree rev_num_always)      # make mpiexec_link when flow123d executable is built
  # explicitly say that source file depends on rev_num_always in order to have actual revision number
  set_source_files_properties(src/main.cc OBJECT_DEPENDS ${CONFIG_FILE_REV_NUM}) 
  set_source_files_properties(src/main.cc OBJECT_DEPENDS ${CONFIG_FILE_HEADER}) 
else()
  message(STATUS "PETSC not found. Omitting configuration of the target 'flow123d'.")
endif()

###############################################################################
# Summary

  
message(STATUS "=============== Summary ===================")
message(STATUS "Boost:      ${Boost_LIBRARIES}")
message(STATUS "Petsc:      ${PETSC_LIBRARIES}")
message(STATUS "   + :      ${PETSC_ADDITIONAL_LIBS}")
message(STATUS "BDDCML:     ${BDDCML_LIBRARIES}")
message(STATUS "Armadillo:  ${ARMADILLO_LIBRARIES}")
message(STATUS "Python:     ${PYTHON_INCLUDE_DIRS}")
message(STATUS "YamlCpp:    ${YamlCpp_LIBRARY}")
message(STATUS "PugiXml:    ${PugiXml_LIBRARY}")
message(STATUS "ZLib:       ${Zlib_LIBRARY}")
message(STATUS "===========================================")
message(STATUS "INCLUDE_DIRECTORIES:")
get_property(dirs DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR} PROPERTY INCLUDE_DIRECTORIES)
foreach(dir ${dirs})
  message(STATUS " - ${dir}")
endforeach()
message(STATUS "===========================================")

####################################################################################
# Install and Packaging

# Cache PLATFORM_NAME
# set(CMAKE_INSTALL_MESSAGE ALWAYS)
set(PLATFORM_NAME ${PLATFORM_NAME} CACHE INTERNAL "Platform name used by prepare_to_publish script.")
set(CMAKE_INSTALL_PREFIX "$ENV{FLOW123D_HOME}") # install prefix is taken from sys env variable FLOW123D_HOME
set(_CPACK_INSTALL_ROOT ${CMAKE_BINARY_DIR}/_CPack_Packages/Linux/TGZ/Flow123d-2.0.0_rc-Linux)  


# install bin directory
install(DIRECTORY bin/config                DESTINATION bin/)
install(DIRECTORY bin/ndiff                 DESTINATION bin/)
install(DIRECTORY bin/paraview              DESTINATION bin/)
install(DIRECTORY bin/python                DESTINATION bin/ PATTERN "*.pyc" EXCLUDE)

# binary files
install(PROGRAMS  build_tree/bin/flow123d   DESTINATION bin/)
install(PROGRAMS  bin/runtest               DESTINATION bin/)
install(PROGRAMS  bin/exec_parallel         DESTINATION bin/)
install(PROGRAMS  bin/input_convert.sh      DESTINATION bin/) 
install(PROGRAMS  bin/yaml_converter/yaml_converter.py    DESTINATION bin/yaml_converter)
install(PROGRAMS  bin/yaml_converter/YAMLConverter.py     DESTINATION bin/yaml_converter)
install(PROGRAMS  bin/yaml_converter/change_rules.py      DESTINATION bin/yaml_converter)
install(PROGRAMS  bin/ndiff/ndiff.pl        DESTINATION bin/ndiff)
install(PROGRAMS  build_tree/bin/mpiexec    DESTINATION bin/)

# plain files
install(FILES     bin/README                DESTINATION bin/)
install(FILES     README.md                 DESTINATION .)
install(FILES     CHANGES.md                DESTINATION .)

# install doc
install(FILES     doc/petsc_options_help    DESTINATION doc/)
install(FILES     doc/LICENSE               DESTINATION doc/)
install(FILES     doc/GPL3                  DESTINATION doc/)
install(FILES     doc/README                DESTINATION doc/)
# PDF and htmldoc are marked as OPTIONAL so this version can be installed on a system where 
# latex and pandoc are not available and where container solution could be dangerou (such as HPC resource)
# TODO perhaps by default, PDF and HTML should be required and only if specified in the config.cmake
#      should they be OPTIONAL
install(FILES     ${CMAKE_BINARY_DIR}/doc/reference_manual/flow123d_doc.pdf DESTINATION doc/ OPTIONAL)
install(DIRECTORY ${CMAKE_BINARY_DIR}/htmldoc/html/src/                     DESTINATION doc/htmldoc OPTIONAL)


# install tests dir
find_program (MAKE_EXECUTABLE               NAMES make gmake)
install(CODE "execute_process(COMMAND \"${MAKE_EXECUTABLE}\" clean-tests WORKING_DIRECTORY \"${CMAKE_SOURCE_DIR}\")")


# python libs, sources and config files
install(DIRECTORY src/python/               DESTINATION lib/python/flow123d)
install(FILES     config/host_table.yaml    DESTINATION config/)


# To install MSVC runtime libraries ... ?
INCLUDE(InstallRequiredSystemLibraries)

SET(CPACK_PACKAGE_DESCRIPTION_SUMMARY "Flow123d - subsurface flow and transport simulator")
SET(CPACK_PACKAGE_VENDOR "Technical University of Liberec")
SET(CPACK_PACKAGE_VERSION_MAJOR ${GIT_VERSION_MAJOR})
SET(CPACK_PACKAGE_VERSION_MINOR ${GIT_VERSION_MINOR})
SET(CPACK_PACKAGE_VERSION_PATCH ${GIT_VERSION_PATCH})
SET(CPACK_COMPONENTS_ALL_IN_ONE_PACKAGE On)
SET(CPACK_SOURCE_IGNORE_FILES "/build_tree/" "/\\\\.git/" "/unit_tests/.*/makefile")
SET(CPACK_PACKAGE_DESCRIPTION_FILE "${CMAKE_CURRENT_SOURCE_DIR}/README.md")
set(CPACK_GENERATOR "TGZ")
set(CPACK_SOURCE_GENERATOR "TGZ")

add_custom_target(package_debug
              COMMAND ${CMAKE_CPACK_COMMAND} --debug --verbose --config CPackConfig.cmake
              COMMENTS "Build package with debug and verbose output")

INCLUDE(CPack)

FILE(APPEND ${CMAKE_BINARY_DIR}/CPackConfig.cmake "set(CPACK_WIN_BINARY_DIRECTORY \"${WIN_CMAKE_CURRENT_BINARY_DIR}\")")
