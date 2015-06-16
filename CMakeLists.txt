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
if (EXISTS "${CMAKE_SOURCE_DIR}/${CONFIG_CMAKE}")
    include("${CMAKE_SOURCE_DIR}/${CONFIG_CMAKE}")
else()
    if (EXISTS "${CMAKE_BINARY_DIR}/${CONFIG_CMAKE}")
        include("${CMAKE_BINARY_DIR}/${CONFIG_CMAKE}")
    else()
        # create from template
        MESSAGE(FATAL_ERROR "Need configuration file: ${CONFIG_CMAKE}. See config.cmake.template.")
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


# Project name and language(s)
project(Flow123d CXX C Fortran)

# CMake version required for pre-build this project
cmake_minimum_required(VERSION 2.8.8)
if (${CMAKE_VERSION} VERSION_GREATER 2.8.4)
  cmake_policy(VERSION 2.8.4)
  cmake_policy(SET CMP0017 NEW)
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




if(FLOW_CC_FLAGS) 
  # user given flags overrides everything  
  message(STATUS "Using user specified compiler options: ${FLOW_CC_FLAGS}")
  set(CMAKE_CXX_FLAGS  "${FLOW_CC_FLAGS}")
  set(CMAKE_C_FLAGS  "${FLOW_CC_FLAGS}")
  set(CMAKE_EXE_LINKER_FLAGS "${FLOW_CC_FLAGS}")
  
else(FLOW_CC_FLAGS)
 
  message(STATUS "Using build type: ${FLOW_BUILD_TYPE}")
  if(FLOW_BUILD_TYPE STREQUAL "debug")
    # debug mode
    set(CMAKE_CXX_FLAGS "-g -O0")
    set(CMAKE_C_FLAGS "-g -O0")
    set(CMAKE_EXE_LINKER_FLAGS "-g ")
    flow_define(DEBUG)
  elseif(FLOW_BUILD_TYPE STREQUAL "release")
    # Release mode
    set(CMAKE_CXX_FLAGS "-O3")
    set(CMAKE_C_FLAGS "-O3")
    # armadillo is configured using this definition
    add_definitions(-DARMA_NO_DEBUG)    # no bound checks in armadillo 
    flow_define(DEBUG_PROFILER)
    flow_define(RUN_UNIT_BENCHMARKS)    # build benchmark unit tests
  elseif(FLOW_BUILD_TYPE STREQUAL "release_with_debug")
    # Release mode with debug info 
    set(CMAKE_CXX_FLAGS "-O3")
    set(CMAKE_C_FLAGS "-O3")
    flow_define(DEBUG)
  elseif(FLOW_BUILD_TYPE STREQUAL "coverage")
    # Coverage mode creates instrumented binary
    SET(CMAKE_CXX_FLAGS "-g -O0 -fprofile-arcs -ftest-coverage")
    SET(CMAKE_C_FLAGS "-g -O0 -fprofile-arcs -ftest-coverage")
    SET(CMAKE_EXE_LINKER_FLAGS "-g -O0 -fprofile-arcs -ftest-coverage")

    # set external libs flags
    SET(EXTERNAL_LIBS_CXX_FLAGS "-g -O0")
    SET(EXTERNAL_LIBS_C_FLAGS "-g -O0")
    flow_define(DEBUG)
    INCLUDE(CodeCoverage)
    # add custom target to collect coverage into xml file
    COVERAGE_SETUP_COBERTURA_COLLECTION (cobertura-collect-coverage coverage.xml)
  else()
  endif()
  
endif(FLOW_CC_FLAGS)


# Compiler-specific C++11 activation.
if ("${CMAKE_CXX_COMPILER_ID}" MATCHES "GNU")        
      execute_process(
          COMMAND ${CMAKE_CXX_COMPILER} -dumpversion OUTPUT_VARIABLE GCC_VERSION)
      set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wno-unused-local-typedefs")

      # is cygwin env?
      if (CYGWIN)
          # is version 4.8 or higher?
          if (GCC_VERSION VERSION_GREATER 4.8 OR GCC_VERSION VERSION_EQUAL 4.8)
            set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=gnu++11")
            flow_define(HAVE_CYGWIN)
          # version is too low
          else ()
            message(SEND_ERROR "Error: Minimal GCC version under Cygwin is 4.8!." )
          endif ()

      # not cygwin 
      # is version 4.7 or higher?
      elseif (GCC_VERSION VERSION_GREATER 4.7 OR GCC_VERSION VERSION_EQUAL 4.7)
          set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")
      # version is too low
      else()
          message(SEND_ERROR "Error: Minimal GCC version is 4.7!.") 
      endif()   

    # Necessary for FParser - TODO: move to better place
    # ?? what about Clang
    flow_define(HAVE_SINCOS)
    
elseif ("${CMAKE_CXX_COMPILER_ID}" MATCHES "Clang")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -stdlib=libc++") # for some reason BOOST doesn't work with option -std=c++11
    
    # BOOST must be compiled with clang in order to use C++11 standard. 
    set(Boost_FORCE_REBUILD "true")
    set(Boost_CONFIG_FLAGS " --with-toolset=clang ")
    set(Boost_BUILD_FLAGS " toolset=clang cxxflags=\"-stdlib=libc++\" linkflags=\"-stdlib=libc++\" " )        
    
    # 
    # TODO: make test for presence of sincos() function and set FLOW123D_HAVE_SINCOS appropriatly, used in fparser extsrc/fpaux.cc
    #
else ()
    message(STATUS "Your C++ compiler does not support C++11. Some features may be disabled.")
endif()


#########################################################################
# Determining timer with highest resolution
include(ProfilerTimerCheck)
set(TIMER_RESOLUTION No)

# try to get best timer on windows
if (WIN32 OR CYGWIN)
    # on windows and cygwin try to get access QueryPerformanceCounter timer
    TEST_TIMER (HAVE_TIMER_QUERY_PERFORMANCE_COUNTER
      ${CMAKE_SOURCE_DIR}/tests/timer_queryperformance.cpp) 

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

message(STATUS "\n =======================================================")
set(Boost_USE_STATIC_LIBS       ON)
set(Boost_USE_MULTITHREADED     OFF)
set(Boost_USE_STATIC_RUNTIME    OFF)
set(EXTERNAL_BOOST_DIR  "${EXTERNAL_PROJECT_DIR}/boost_build")  
set(REQUIRED_Boost_COMPONENTS program_options serialization regex filesystem system)
 
set(Boost_DEBUG TRUE)
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

message(STATUS "BOOST_ROOT: ${BOOST_ROOT}")
message(STATUS "Boost_LIBRARIES: ${Boost_LIBRARIES}")
message(STATUS "Boost_LIBRARY_DIRS: ${Boost_LIBRARY_DIRS}")   
message(STATUS "Boost_INCLUDE_DIR: ${Boost_INCLUDE_DIR}")


# We are looking for PETSC, f succed this module set:
#  PETSC_FOUND        - system has PETSc
#  PETSC_INCLUDES     - the PETSc include directories
#  PETSC_LIBRARIES    - Link these to use PETSc
#  PETSC_COMPILER     - Compiler used by PETSc, helpful to find a compatible MPI
#  PETSC_DEFINITIONS  - Compiler switches for using PETSc
#  PETSC_MPIEXEC      - Executable for running MPI programs
#  PETSC_VERSION      - Version string (MAJOR.MINOR.SUBMINOR)

message(STATUS "\n =======================================================")
include (FindPackageMultipass)

# this export list is for possible BDDCML installation
set(PETSC_EXPORT_LIST 
  MPIFC_SHOW MPICXX_SHOW MPICC_SHOW CC FC CC_FLAGS PCC_FLAGS FC_FLAGS CPP_FLAGS
  METIS_INCLUDE BLASLAPACK_INCLUDE PARMETIS_INCLUDE MUMPS_INCLUDE
  MPI_LIB MUMPS_LIB  METIS_LIB PARMETIS_LIB BLASLAPACK_LIB SCALAPACK_LIB
  UMFPACK_LIB)  
find_package(PETSc 3.4.0)


if (NOT PETSC_FOUND)
    if (PETSC_DIR) 
        message(STATUS "Can not configure with given PETSC_DIR=${PETSC_DIR}")
    else()    
        # Install PETSC
        # also set PETSC_DIR, PETSC_ARCH
        set(PETSC_INSTALL_DIR  ${EXTERNAL_PROJECT_DIR}/petsc_build)
        set(PETSC_INSTALL_URL \"http://bacula.nti.tul.cz/~jan.brezina/flow123d_libraries/petsc-lite-3.4.3.tar.gz\")
        include(InstallPetsc)
        
        set(PETSC_CURRENT "NO")
        find_package(PETSc 3.4.0)
        
    endif()    
        
       
    # If we haven't PETSC, we try to find MPI
        
    # === Variables ===
    # This module will set the following variables per language in your project,
    # where <lang> is one of C, CXX, or Fortran:
    #   MPI_<lang>_FOUND           TRUE if FindMPI found MPI flags for <lang>
    #   MPI_<lang>_COMPILER        MPI Compiler wrapper for <lang>
    #   MPI_<lang>_COMPILE_FLAGS   Compilation flags for MPI programs
    #   MPI_<lang>_INCLUDE_PATH    Include path(s) for MPI header
    #   MPI_<lang>_LINK_FLAGS      Linking flags for MPI programs
    #   MPI_<lang>_LIBRARIES       All libraries to link MPI programs against
    # Additionally, FindMPI sets the following variables for running MPI
    # programs from the command line:
    #   MPIEXEC                    Executable for running MPI programs
    #   MPIEXEC_NUMPROC_FLAG       Flag to pass to MPIEXEC before giving
    #                              it the number of processors to run on
    #   MPIEXEC_PREFLAGS           Flags to pass to MPIEXEC directly
    #                              before the executable to run.
    #   MPIEXEC_POSTFLAGS          Flags to pass to MPIEXEC after other flags
#
#    find_package(MPI)
#  
#    if (MPI_CXX_FOUND)
#        flow_define(HAVE_MPI)
#    else(MPI_CXX_FOUND) 
        # build single processor replacement taken from PETSc
#        set(mpiuni_sources    system/mpiuni/mpi.c)
#        set(MPIUNI_INCLUDES src/system/mpiuni)
#        flow_define(HAVE_MPI)
#        flow_define(HAVE_MPIUNI)
#    
        # overwrite NOTFOUND status
#        set(MPI_CXX_LIBRARIES mpiuni_lib)
#        set(MPI_CXX_INCLUDE_PATH ${MPIUNI_INCLUDES})
#    
#    endif(MPI_CXX_FOUND)

endif(NOT PETSC_FOUND)

if(PETSC_FOUND)
    flow_define(HAVE_PETSC)
    flow_define(HAVE_MPI)
    set(MPI_CXX_LIBRARIES ${PETSC_EXTERNAL_LIB})
    set(MPI_CXX_INCLUDE_PATH ${PETSC_INCLUDES})
    
    message(STATUS "PETSC_LIBS: ${PETSC_LIBRARIES}")
    message(STATUS "PETSC_EXTERNAL_LIB: ${PETSC_EXTERNAL_LIB}")
    message(STATUS "PETSC_EXTERNAL_LIB: ${petsc_libraries_external}")
endif(PETSC_FOUND)

#find_package(PARMETIS)

message(STATUS "\n =======================================================")
#################################################################################
#  BDDCML_FOUND - set to true if the library is found

if (BDDCML_ROOT) 
    find_package(BDDCML)

    if (NOT BDDCML_FOUND)
        set(EXTERNAL_BLOPEX_DIR  ${EXTERNAL_PROJECT_DIR}/blopex_build)
        include(InstallBLOPEX)
        set(EXTERNAL_BDDCML_DIR  ${EXTERNAL_PROJECT_DIR}/bddcml_build)
        include(InstallBDDCML)

        set(BDDCML_bddcml_config_current "NO")
        find_package(BDDCML)
    endif()
    
    if (BDDCML_FOUND)
        flow_define(HAVE_BDDCML)
    endif()
else()
    # set(BDDCML_LIBRARIES "")
endif()    


#################################################################################
#  Armadillo_FOUND - set to true if the library is found
#  Armadillo_INCLUDE_DIRS - list of required include directories
#  Armadillo_LIBRARIES - list of libraries to be linked 
#  Armadillo_CONFIGURE_CMAKE - path to configure cmake script that introduce target 'armadillo'
#                              with correctly set dependencies on other libraries

message(STATUS "\n =======================================================")
find_package(ArmadilloFlow 3.2.4)
if (NOT Armadillo_FOUND)
  message(STATUS "Try to install our own Armadillo.")
  set(EXTERNAL_ARMADILLO_DIR  ${EXTERNAL_PROJECT_DIR}/armadillo_build)
  include(InstallArmadillo)

  # we try to erease from the cache everything that can cause failure of the second find_package
  # however just unsetting Boost_INCLUDE_DIR should be enough
  unset(Armadillo_FOUND CACHE) 
  unset(Armadillo_INCLUDE_DIR CACHE)        # overwrite cache    
  unset(Armadillo_LIBRARY_DIRS CACHE)
  unset(Armadillo_CONFIGURE_FILE CACHE)

  set(Armadillo_ROOT_HINT ${EXTERNAL_ARMADILLO_DIR})  
  find_package(ArmadilloFlow 3.2.4 REQUIRED)

endif()


message(STATUS "Armadillo_LIBRARIES: ${Armadillo_LIBRARIES}")
message(STATUS "Armadillo_INCLUDE_DIRS: ${Armadillo_INCLUDE_DIRS}")
message(STATUS "Armadillo_CONFIGURE_FILE: ${Armadillo_CONFIGURE_FILE}")
include("${Armadillo_CONFIGURE_FILE}")   # this introduce target "armadillo" with correct dependencies
IF (PETSC_FOUND)
    # until we know how to take only blas and lapack from PETSC we have to make Armadillo depend on whole PETSC lib.
    # this may be fragile because of _NOCONFIG
    set_target_properties(armadillo PROPERTIES 
      IMPORTED_LINK_INTERFACE_LIBRARIES_NOCONFIG "${PETSC_EXTERNAL_LIB}"
      ) 
endif()

####################################################################################
# PYTHON
#
# optionaly we use Python interface
if (USE_PYTHON)
    message(STATUS "\n ======================================================= ")
    find_package(PythonLibs 2.7)
    if (PYTHONLIBS_FOUND)
        include_directories( ${PYTHON_INCLUDE_DIRS})
        flow_define(HAVE_PYTHON)
        
        string(REGEX REPLACE "([0-9]+).*" "\\1" PYTHONLIBS_VERSION_MAJOR ${PYTHONLIBS_VERSION_STRING})
        string(REGEX REPLACE "[0-9]+.([0-9]+).*" "\\1" PYTHONLIBS_VERSION_MINOR ${PYTHONLIBS_VERSION_STRING})
        flow_define_constant(PYTHONLIBS_VERSION_MAJOR ${PYTHONLIBS_VERSION_MAJOR})
        flow_define_constant(PYTHONLIBS_VERSION_MINOR ${PYTHONLIBS_VERSION_MINOR})
        flow_define(PYTHONLIBS_VERSION_STRING ${PYTHONLIBS_VERSION_STRING})
        message(STATUS "Python version: ${PYTHONLIBS_VERSION_STRING}\n")
        message(STATUS "Python major version: ${PYTHONLIBS_VERSION_MAJOR}\n")
        message(STATUS "Python minor version: ${PYTHONLIBS_VERSION_MINOR}\n")
        message(STATUS "Python include: ${PYTHON_INCLUDE_DIRS}")
        message(STATUS "Python libs: ${PYTHON_LIBRARIES}")

        # set extra python module path in src directory, if case PYTHON_COPY if OFF
        set(PYTHON_EXTRA_MODULES_PATH "${CMAKE_SOURCE_DIR}/src/python")
        
        if (PYTHON_COPY)
            set(TARGET_PREFIX "${CMAKE_BINARY_DIR}")
            set(PYTHON_EXTRA_MODULES_PATH "${TARGET_PREFIX}/lib/flow123d")

            string(REGEX MATCH "python[^/]*$" PYTHON_SUBDIR ${PYTHON_INCLUDE_DIRS})
            string(REGEX MATCH "python[0-9]\\.[0-9]" PYTHON_SUBDIR_MINOR ${PYTHON_SUBDIR})
            
            # copy python files to build tree 
            # first get prefix (we support just one prefix)
            execute_process(
                COMMAND python -c "import sys; print sys.prefix"
                OUTPUT_VARIABLE PY_PREFIX
                OUTPUT_STRIP_TRAILING_WHITESPACE)        
            execute_process(
                COMMAND python -c "import sys; print sys.exec_prefix"
                OUTPUT_VARIABLE PY_EXEC_PREFIX
                OUTPUT_STRIP_TRAILING_WHITESPACE)        
            if (NOT (${PY_PREFIX} STREQUAL ${PY_EXEC_PREFIX}) )
                message(WARNING "Python exec_prefix (${PY_EXEC_PREFIX}) different then prefix (${PY_PREFIX}). Not supported.")
            endif()
            message(STATUS "Python prefix: ${PY_PREFIX}")
            message(STATUS "Python target prefix: ${TARGET_PREFIX}")
            flow_define(PYTHON_PREFIX ${PY_PREFIX})
            
            # copy python library
            set(PYTHON_COPY_LIB_TARGET  lib/${PYTHON_SUBDIR_MINOR})
            execute_process(
                COMMAND mkdir -p ${TARGET_PREFIX}/lib 
                )             
            # recursive copy, not overwrite existing files
            # we copy all and exclude some paths and files during install    
            execute_process(
                COMMAND cp -rn ${PY_PREFIX}/lib/${PYTHON_SUBDIR_MINOR} ${TARGET_PREFIX}/lib 
                )
            # copy flow123d python scripts and modules
            execute_process(
                COMMAND cp -rT ${CMAKE_SOURCE_DIR}/src/python ${TARGET_PREFIX}/lib/flow123d 
                )

            # copy python include            
            set(PYTHON_COPY_INCLUDE_TARGET include/${PYTHON_SUBDIR})
            execute_process(          
                COMMAND mkdir -p ${TARGET_PREFIX}/include 
                )             
            # recursive copy, not overwrite existing files
            # we copy all and exclude some paths and files during install    
            execute_process(
                COMMAND cp -rn ${PYTHON_INCLUDE_DIR} ${TARGET_PREFIX}/include 
                )                           
            
                
        endif()
    else()
        set(PYTHON_INCLUDE_DIRS "")
        set(PYTHON_LIBRARIES "")
    endif()

    # define varible for storing this extra module path
    # this variable is used later to update search module path
    # if PYTHON_COPY is present it points to lib/flow123d dir in build tree
    # if not, default location is src/python dir in CMake_SOURCE_DIR location
    flow_define(PYTHON_EXTRA_MODULES_PATH "${PYTHON_EXTRA_MODULES_PATH}")
endif()  

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
    ${BDDCML_INCLUDES}
    ${Boost_INCLUDE_DIR}
    ${Armadillo_INCLUDE_DIRS}
    ${MPI_CXX_INCLUDE_PATH}
    ${CMAKE_BINARY_DIR}/src/dealii/include
    ${CMAKE_SOURCE_DIR}/src/dealii/include
    ${CMAKE_SOURCE_DIR}/third_party/tbb43_20150316oss/include
)

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
      message(STATUS "using petsc_mpiexec: ${RESOLVED_MPIEXEC}")
      set(MPIEXEC_PATH "${RESOLVED_MPIEXEC}")
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

configure_file(${CMAKE_SOURCE_DIR}/CMake/mpiexec_link_template ${CMAKE_BINARY_DIR}/bin/mpiexec)


##################################################################################################################### 
# Macro to make usefull links in SOURCE_DIR
##################################


# initialize a global variable (cached)
set(SOURCE_TREE_LINKS  "" CACHE INTERNAL "" FORCE)
              

macro(make_link_in_source FILE)
  # get relative filename if absolute address is given
  string(REPLACE "${CMAKE_SOURCE_DIR}/" "" rel_path "${FILE}")
  execute_process(
    COMMAND ${CMAKE_COMMAND} -E create_symlink ${CMAKE_SOURCE_DIR}/build_tree/${rel_path} ${CMAKE_SOURCE_DIR}/${rel_path}
    WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
  )
  
  set(SOURCE_TREE_LINKS  ${SOURCE_TREE_LINKS} ${CMAKE_SOURCE_DIR}/${FILE} CACHE INTERNAL "")
endmacro()

make_link_in_source(bin/flow123d)
make_link_in_source(bin/mpiexec)
make_link_in_source(bin/tests/Makefile)
make_link_in_source(unit_tests/Makefile)

############################################################33 process subdirectories


# add fparser sources
message(STATUS "add subdir FParser")
add_subdirectory(third_party/fparser-4.5.1)

# add tbb sources
message(STATUS "add subdir TBB")
add_subdirectory(third_party/tbb43_20150316oss)

# add various our sources
add_subdirectory(src)
add_subdirectory(doc)
add_subdirectory(unit_tests) # this also adds third_party/gtest

# optionally make targets for building reference_manual
# we have to call separate instance of CMake since UseLatex.cmake
# depends on many tools. If any of them miss CMake fails, while we can live
# without targets for Latex doc.
message(STATUS "==== child CMake for Latex sources ====") 
execute_process(COMMAND ${CMAKE_COMMAND} -E make_directory ${CMAKE_BINARY_DIR}/doc/reference_manual)
execute_process(COMMAND ${CMAKE_COMMAND} ${CMAKE_SOURCE_DIR}/doc/reference_manual
                WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/doc/reference_manual
                RESULT_VARIABLE res_var)
message(STATUS "==== back in main CMake ====")
if(NOT "${res_var}" STREQUAL "0")
  message(WARNING "Some problem with Latex, 'make ref-doc' may not work.")
endif()


#################################################################### add targets


if( PETSC_FOUND )
  ########################
# Setting for RPATH

 
  SET (CMAKE_EXE_LINKER_FLAGS
    "${CMAKE_EXE_LINKER_FLAGS} -Wl,-rpath,$ORIGIN/../lib" )
  
  
  ### add target for flow123d binary
  add_executable(bin/flow123d src/main.cc)
  target_link_libraries(bin/flow123d  flow123d_lib  ${PETSC_ADDITIONAL_LIBS}  semchem  ${PYTHON_LIBRARIES} ${Boost_LIBRARIES})     # set libraries linked to main binary (order matter)
  add_dependencies(bin/flow123d update_source_tree rev_num_always mpiexec_link)      # make mpiexec_link when flow123d executable is built
  # explicitly say that source file depends on rev_num_always in order to have actual revision number
  set_source_files_properties(src/main.cc OBJECT_DEPENDS ${CONFIG_FILE_REV_NUM}) 
  set_source_files_properties(src/main.cc OBJECT_DEPENDS ${CONFIG_FILE_HEADER}) 
else()
  message(STATUS "PETSC not found. Omitting configuration of the target 'flow123d'.")
endif()

###############################################################################
# Summary

  
message(STATUS "=============== Summary ===================")
message(STATUS "Boost:      ${BOOST_LIBRARIES}")
message(STATUS "Petsc:      ${PETSC_LIBRARIES}")
message(STATUS "   + :      ${PETSC_ADDITIONAL_LIBS}")
message(STATUS "BDDCML:     ${BDDCML_LIBRARIES}")
message(STATUS "Armadillo:  ${ARMADILLO_LIBRARIES}")
message(STATUS "Python:     ${PYTHON_INCLUDE_DIRS}")
message(STATUS "===========================================")


#include(EchoTargets)
#message(STATUS "Lib sufix: ${CMAKE_FIND_LIBRARY_SUFFIX}")
#echo_targets("bin/flow123d")

####################################################################################
# Install and Packaging

# Cache PLATFORM_NAME
set(PLATFORM_NAME ${PLATFORM_NAME} CACHE INTERNAL "Platform name used by prepare_to_publish script.")

if (UNIX)
  set(CMAKE_INSTALL_PREFIX "$ENV{HOME}/local/flow123d")
endif()


# install bin directory
install(DIRECTORY bin/config              DESTINATION bin/)
install(DIRECTORY bin/ndiff               DESTINATION bin/)
install(DIRECTORY bin/paraview            DESTINATION bin/)
install(DIRECTORY bin/python                  DESTINATION bin/)
install(PROGRAMS  bin/tests/run_test.sh   DESTINATION bin/tests)

install(PROGRAMS  bin/flow123d.sh         DESTINATION bin/)
install(PROGRAMS  bin/ini2json.sh         DESTINATION bin/)
install(PROGRAMS  bin/time_limit.sh       DESTINATION bin/)

# install shared libraries
if (CYGWIN)
    install(CODE "execute_process(COMMAND \"${CMAKE_SOURCE_DIR}/bin/copy_shared_libs.sh\" cygwin WORKING_DIRECTORY \"${CMAKE_BINARY_DIR}\")")
    install(DIRECTORY build_tree/bin/            
        DESTINATION bin/
        PATTERN mpiexec EXCLUDE)

else()
    install(CODE "execute_process(COMMAND \"${CMAKE_SOURCE_DIR}/bin/copy_shared_libs.sh\" linux WORKING_DIRECTORY \"${CMAKE_BINARY_DIR}\")")
    install(DIRECTORY build_tree/lib/            
        DESTINATION lib/
        PATTERN lib/python* EXCLUDE
        PATTERN *.a EXCLUDE
        )
endif()        

# install python modules
if (PYTHON_COPY)
    install(DIRECTORY build_tree/${PYTHON_COPY_LIB_TARGET} 
            DESTINATION lib
            PATTERN *.dll EXCLUDE
            PATTERN *.a EXCLUDE
            PATTERN *.so EXCLUDE
            PATTERN */dist-packages EXCLUDE
            PATTERN *.py EXCLUDE)
    install(DIRECTORY build_tree/${PYTHON_COPY_INCLUDE_TARGET} 
            DESTINATION include)
endif()

install(TARGETS    bin/flow123d  
  RUNTIME DESTINATION bin)
install(PROGRAMS   ${MPIEXEC_PATH}              DESTINATION bin/)
install(FILES bin/README                        DESTINATION bin/)

# make manual optional until we have Latex on Windows build machine
install(DIRECTORY doc/con_syntax                DESTINATION doc/)
install(FILES ${CMAKE_BINARY_DIR}/doc/reference_manual/flow123d_doc.pdf 
                                                DESTINATION doc/)
install(FILES doc/petsc_options_help            DESTINATION doc/)
install(FILES doc/LICENSE                       DESTINATION doc/)
install(FILES doc/GPL3                          DESTINATION doc/)
install(FILES doc/README                        DESTINATION doc/)


# clean tests directory before install
find_program (MAKE_EXECUTABLE NAMES make gmake)
install(CODE "execute_process(COMMAND \"${MAKE_EXECUTABLE}\" clean-tests WORKING_DIRECTORY \"${CMAKE_SOURCE_DIR}\")")
install(DIRECTORY tests DESTINATION .)


# copy base directory
install(FILES README.md DESTINATION .)
install(FILES CHANGES.md DESTINATION .)

# To install MSVC runtime libraries ... ?
INCLUDE(InstallRequiredSystemLibraries)


SET(CPACK_PACKAGE_DESCRIPTION_SUMMARY "Flow123d - subsurface flow and transport simulator")
SET(CPACK_PACKAGE_VENDOR "Technical University of Liberec")
SET(CPACK_PACKAGE_VERSION_MAJOR ${GIT_VERSION_MAJOR})
SET(CPACK_PACKAGE_VERSION_MINOR ${GIT_VERSION_MINOR})
SET(CPACK_PACKAGE_VERSION_PATCH ${GIT_VERSION_PATCH})

SET(CPACK_SOURCE_IGNORE_FILES "/build_tree/" "/\\\\.git/" "/unit_tests/.*/makefile")

macro(MAKE_WIN_PATH VAR)
    execute_process(
        COMMAND cygpath -aw "${${VAR}}"
        OUTPUT_VARIABLE CYGPATH_${VAR}
        OUTPUT_STRIP_TRAILING_WHITESPACE)        
    string(REGEX REPLACE "\\\\" "\\\\\\\\" WIN_${VAR} ${CYGPATH_${VAR}}) 
endmacro()

if (CYGWIN) 
    ## Try to use NSIS packager ##
    #
    # Package production runs as follows:
    # 1) running CMake the file CPackConfig.cmake is created and a target 'package'
    # 2) make package runs CPack which creates NSIS script file based on template NSIS.template.in in the Module path
    # 3) CPack run NSIS to compress package and create the installer
    #
    # Following hacks are necessary to combine CYGWIN and NSIS packager.
    # 1) We add modified NSIS.template.in into our modules directory CMake/Modules and 
    #    set this as the first one CMAKE_MODULE_PATH. 
    # 2) We have to convert POSIX path to binary directory and use this 
    #    to set all paths for NSIS installer
    # 3) We pass converted BINARY_DIR path to CPack, adding it manually to the end of 
    #    CPackConfig.cmake generated by CMake (must proceed after include(CPack))
    # 4) The passed variable is used in our NSIS.template.in to set working directory and directory of the result.
    # 5) CPack call makensis.exe with POSIX path to the script. So this end with error. We add target 'post_package' that calls NSIS with converted path.
    #    This target has to by run manually.
    #    
    # We need Windows style paths.
    MAKE_WIN_PATH(CMAKE_CURRENT_SOURCE_DIR)
    MAKE_WIN_PATH(CMAKE_CURRENT_BINARY_DIR)      
    
    
    # 1. Page
    # Displayed name of the software, and target subdirectory of installation prefix.
    SET(CPACK_PACKAGE_INSTALL_DIRECTORY "Flow123d ${GIT_VERSION_MAJOR}.${GIT_VERSION_MINOR}.${GIT_VERSION_PATCH}")
    
    # 2. Page
    # Software banner in installer.
    SET(CPACK_PACKAGE_ICON "${WIN_CMAKE_CURRENT_SOURCE_DIR}\\\\bin\\\\graphics\\\\flow123d.ico")
    # Displayed license  agreement.
    SET(CPACK_RESOURCE_FILE_LICENSE "${WIN_CMAKE_CURRENT_SOURCE_DIR}\\\\doc\\\\LICENSE")
    SET(CPACK_PACKAGE_DESCRIPTION_FILE "${WIN_CMAKE_CURRENT_SOURCE_DIR}\\\\README.md")
    
    # 3. Page
    # Modify PATH toggle. If this is set to "ON", then an extra page
    # will appear in the installer that will allow the user to choose
    # whether the program directory should be added to the system PATH variable.
    set(CPACK_NSIS_MODIFY_PATH "ON")
    set(CPACK_CREATE_DESKTOP_LINKS flow123d)
    set(CPACK_PACKAGE_EXECUTABLES "flow123d" "Flow123d")
    # Specify links in [application] menu. This should contain a list of pair 
    # "link" "link name". The link may be an URL or a path relative to installation 
    # prefix. Like:
    #               set(CPACK_NSIS_MENU_LINKS
    #                   "doc/cmake-@CMake_VERSION_MAJOR@.@CMake_VERSION_MINOR@/cmake.html" "CMake Help"
    #                   "http://www.cmake.org" "CMake Web Site")
    set(CPACK_NSIS_MENU_LINKS 
        "flow123d" "Flow123d"
        "http://flow123d.github.io" "Flow123d web page")
    
    
    
    # 4. Page
    # The default installation directory presented
    # to the end user by the NSIS installer is under this root dir. The full
    # directory presented to the end user is:
    # ${CPACK_NSIS_INSTALL_ROOT}/${CPACK_PACKAGE_INSTALL_DIRECTORY}
    #set(CPACK_NSIS_INSTALL_ROOT "")
    
    
    
    # An icon filename. The name of a *.ico file used as the main icon for the generated
    # install program.
    #set(CPACK_NSIS_MUI_ICON "")

    # An icon filename. The name of a *.ico file used as the main icon for the generated
    # uninstall program.
    #set(CPACK_NSIS_MUI_UNIICON "")

    # Extra NSIS commands that will be added to the beginning of the install Section,
    # before your install tree is available on the target system.
    #set(CPACK_NSIS_EXTRA_PREINSTALL_COMMANDS "")

    # Extra NSIS commands that will be added to the end of the install Section, 
    # after your install tree is available on the target system.
    
    # File association for CON files."
    set(CPACK_NSIS_EXTRA_INSTALL_COMMANDS "
        WriteRegStr HKCR '.con' '' 'Flow123dFile'
        WriteRegStr HKCR 'Flow123dFile\\\\shell' '' 'open'
        WriteRegStr HKCR 'Flow123dFile\\\\shell\\\\open\\\\command' '' '$INSTDIR\\\\bin\\\\flow123d.exe \\"%1\\"'
        System::Call 'Shell32::SHChangeNotify(i 0x8000000, i 0, i 0, i 0)'
        ")

    # Extra NSIS commands that will be added to the uninstall Section,
    # before your install tree is removed from the target system.
    set(CPACK_NSIS_EXTRA_UNINSTALL_COMMANDS "
        DeleteRegKey HKCR '.con'
        DeleteRegKey HKCR 'Flow123dFile'
        ")

    # The arguments that will be passed to the NSIS SetCompressor command.
    #set(CPACK_NSIS_COMPRESSOR "")

    # Ask about uninstalling previous versions first.
    # If this is set to "ON", then an installer will look for previous
    # installed versions and if one is found, ask the user whether to
    # uninstall it before proceeding with the install.
    #set(CPACK_NSIS_ENABLE_UNINSTALL_BEFORE_INSTALL "")

    
    # A path to the executable that contains the installer icon.
    #set(CPACK_NSIS_INSTALLED_ICON_NAME "")

    # URL to a web site providing assistance in installing your application.
    set(CPACK_NSIS_HELP_LINK "http://flow123d.github.io")

    # URL to a web site providing more information about your application.
    set(CPACK_NSIS_URL_INFO_ABOUT "http://flow123d.github.io")

    # Contact information for questions and comments about the installation process.
    set(CPACK_NSIS_CONTACT "jan.brezin@tul.cz")

    # Additional NSIS commands for creating start menu shortcuts.
    set(CPACK_NSIS_CREATE_ICONS_EXTRA "")

    # Additional NSIS commands to uninstall start menu shortcuts.
    #set(CPACK_NSIS_DELETE_ICONS_EXTRA "")

    # Creating NSIS start menu links assumes that they are in 'bin' unless this 
    # variable is set. For example, you would set this to 'exec' if your executables 
    # are in an exec directory.
    #set(CPACK_NSIS_EXECUTABLES_DIRECTORY "")

    # Specify an executable to add an option to run on the finish page 
    # of the NSIS installer.
    #set(CPACK_NSIS_MUI_FINISHPAGE_RUN "")
    

    #SET(CPACK_PACKAGE_ICON "${CMake_SOURCE_DIR}\\\\bin\\\\graphics\\\\flow123d.ico")
    #SET(CPACK_NSIS_DISPLAY_NAME "${CPACK_PACKAGE_INSTALL_DIRECTORY} Flow123d")


    set(CPACK_GENERATOR "NSIS")
    set(CPACK_SOURCE_GENERATOR "ZIP")
    
    add_custom_target(post_package
                      COMMAND makensis.exe "${WIN_CMAKE_CURRENT_BINARY_DIR}\\\\_CPack_Packages\\\\CYGWIN\\\\NSIS\\\\project.nsi"
                      WORKING_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/_CPack_Packages/CYGWIN/NSIS"
                      )
elseif (UNIX)
    SET(CPACK_PACKAGE_DESCRIPTION_FILE "${CMAKE_CURRENT_SOURCE_DIR}/README.md")
    set(CPACK_GENERATOR "TGZ")
    set(CPACK_SOURCE_GENERATOR "TGZ")
    add_custom_target(post_package
                      COMMAND echo ""
                      )
endif()



#IF(WIN32 AND NOT UNIX)
  # There is a bug in NSI that does not handle full unix paths properly. Make
  # sure there is at least one set of four (4) backlasshes.
#  SET(CPACK_PACKAGE_ICON "${CMake_SOURCE_DIR}/Utilities/Release\\\\InstallIcon.bmp")
#  SET(CPACK_NSIS_INSTALLED_ICON_NAME "bin\\\\MyExecutable.exe")
#  SET(CPACK_NSIS_DISPLAY_NAME "${CPACK_PACKAGE_INSTALL_DIRECTORY} My Famous Project")
#  SET(CPACK_NSIS_HELP_LINK "http:\\\\\\\\www.my-project-home-page.org")
#  SET(CPACK_NSIS_URL_INFO_ABOUT "http:\\\\\\\\www.my-personal-home-page.com")
#  SET(CPACK_NSIS_CONTACT "me@my-personal-home-page.com")
#  SET(CPACK_NSIS_MODIFY_PATH ON)
#ELSE(WIN32 AND NOT UNIX)
#  SET(CPACK_STRIP_FILES "bin/MyExecutable")
#  SET(CPACK_SOURCE_STRIP_FILES "")
#ENDIF(WIN32 AND NOT UNIX)
#SET(CPACK_PACKAGE_EXECUTABLES "MyExecutable" "My Executable")

INCLUDE(CPack)

FILE(APPEND ${CMAKE_BINARY_DIR}/CPackConfig.cmake "set(CPACK_WIN_BINARY_DIRECTORY \"${WIN_CMAKE_CURRENT_BINARY_DIR}\")")
