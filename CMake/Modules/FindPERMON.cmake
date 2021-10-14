# - Try to find PERMON
#
#
# Usage:
#  find_package(PERMON)
#
# Setting these changes the behavior of the search
#  PERMON_DIR - directory in which PERMON resides
#  PERMON_EXPORT_LIST - List of variables from PERMON makefile system to export. For each name XYZ in the list
#                      we set CMake variable PERMON_VAR_XYZ.
#  PERMON_ADDITIONAL_LIBS - Add libraries to the link sequence; use only to fix possible reported 
#                          problems with absolute path resolution.
#
# Once done this will define:
#
#  PERMON_FOUND        - system has PERMON
#  PERMON_INCLUDES     - the PERMON include directories
#  PERMON_LIBRARIES    - Link these to use PERMON
#  PERMON_COMPILER     - Compiler used by PERMON, helpful to find a compatible MPI
#  PERMON_DEFINITIONS  - Compiler switches for using PERMON
#  PERMON_MPIEXEC      - Executable for running MPI programs
#  PERMON_VERSION      - Version string (MAJOR.MINOR.SUBMINOR)
#
#  PERMON_EXTERNAL_LIB   - CMake list of resolved (hopefully) external libraries linked by PERMON, 
#                        in the case of static PERMON libraries this list is already included in PERMON_LIBRARIES
#
#  PERMON_VAR_XYZ       - exported variables from given PERMON_EXPORT_LIST

#
# Redistribution and use is allowed according to the terms of the BSD license.
# For details see the accompanying COPYING-CMAKE-SCRIPTS file.
#

if(NOT PERMON_FIND_COMPONENTS)
  set(PERMON_LANGUAGE_BINDINGS "C")
else()
  # Right now, this is designed for compatability with the --with-clanguage option, so
  # only allow one item in the components list.
  list(LENGTH ${PERMON_FIND_COMPONENTS} components_length)
  if(${components_length} GREATER 1)
    message(FATAL_ERROR "Only one component for PERMON is allowed to be specified")
  endif()
  # This is a stub for allowing multiple components should that time ever come. Perhaps
  # to also test Fortran bindings?
  foreach(component ${PERMON_FIND_COMPONENTS})
    list(FIND PERMON_VALID_COMPONENTS ${component} component_location)
    if(${component_location} EQUAL -1)
      message(FATAL_ERROR "\"${component}\" is not a valid PERMON component.")
    else()
      list(APPEND PERMON_LANGUAGE_BINDINGS ${component})
    endif()
  endforeach()
endif()

function (permon_get_version)
  #git describe --tags
  #if (EXISTS  "${PERMON_DIR}")
  #  execute_process(COMMAND  which grpc_cpp_plugin OUTPUT_VARIABLE PERMON_VERSION_MAJOR)
  #  file (STRINGS "${PERMON_DIR}/include/permonversion.h" vstrings REGEX "#define PERMON_VERSION_(RELEASE|MAJOR|MINOR|SUBMINOR|PATCH) ")
  #  foreach (line ${vstrings})
  #    string (REGEX REPLACE " +" ";" fields ${line}) # break line into three fields (the first is always "#define")
  #    list (GET fields 1 var)
  #    list (GET fields 2 val)
  #    set (${var} ${val} PARENT_SCOPE)
  #    set (${var} ${val})         # Also in local scope so we have access below
  #  endforeach ()
  #  if (PERMON_VERSION_RELEASE)
  #    set (PERMON_VERSION "${PERMON_VERSION_MAJOR}.${PERMON_VERSION_MINOR}.${PERMON_VERSION_SUBMINOR}p${PERMON_VERSION_PATCH}" PARENT_SCOPE)
  #  else ()
  #    # make dev version compare higher than any patch level of a released version
  #    set (PERMON_VERSION "${PERMON_VERSION_MAJOR}.${PERMON_VERSION_MINOR}.${PERMON_VERSION_SUBMINOR}.99" PARENT_SCOPE)
  #  endif ()
  #else ()
  #  message (SEND_ERROR "PERMON_DIR can not be used, ${PERMON_DIR} does not exist")
  #endif ()

  #MESSAGE(STATUS "PERMON version: ${PERMON_VERSION}")
  
endfunction ()



########################################################## 
# Try to find main header file.
##########################################################
find_path (PERMON_DIR include/permonsys.h
  HINTS ENV PERMON_DIR
  PATHS
  $ENV{HOME}/permon
  DOC "PERMON Directory")

find_program (MAKE_EXECUTABLE NAMES make gmake)


########################################################## 
# Try to find PERMON_ARCH subdirectory that contains permonconf.h file.
##########################################################
#if (PERMON_DIR AND NOT PERMON_ARCH)
#  set (_permon_arches
#    $ENV{PERMON_ARCH}                   # If set, use environment variable first
#    linux-gnu-c-debug linux-gnu-c-opt  # Debian defaults
#    x86_64-unknown-linux-gnu i386-unknown-linux-gnu)
#  set (permonconf "NOTFOUND" CACHE FILEPATH "Cleared" FORCE)
#  foreach (arch ${_permon_arches})
#    if (NOT PERMON_ARCH)
#      find_path (permonconf permonconf.h
#        HINTS ${PERMON_DIR}
#        PATH_SUFFIXES ${arch}/include bmake/${arch}
#        NO_DEFAULT_PATH)
#      if (permonconf)
#        set (PERMON_ARCH "${arch}" CACHE STRING "PERMON build architecture")
#      endif (permonconf)
#    endif (NOT PERMON_ARCH)
#  endforeach (arch)
#  set (permonconf "NOTFOUND" CACHE INTERNAL "Scratch variable" FORCE)
#endif (PERMON_DIR AND NOT PERMON_ARCH)


########################################################## 
# Use package multipass to alow multiple tries when search for PERMON.
##########################################################
set (permon_slaves LIBRARIES_SYS LIBRARIES_VEC LIBRARIES_MAT LIBRARIES_KSP LIBRARIES_PC
  LIBRARIES_QP LIBRARIES_QPC LIBRARIES_QPPF LIBRARIES_QPS
  INCLUDE_DIR INCLUDE_CONF) #TODO ?
include (FindPackageMultipass)
find_package_multipass (PERMON permon_config_current
  STATES DIR ARCH
  DEPENDENTS INCLUDES LIBRARIES COMPILER MPIEXEC EXTERNAL_LIB ${permon_slaves})
  
#message(STATUS "permon_external_lib: ${PERMON_EXTERNAL_LIB}\n")
#message(STATUS "permon_includes: ${PERMON_INCLUDES}\n")
#message(STATUS "permon_config_current: ${permon_config_current}\n")

  
########################################################## 
# Detect 'rules' and 'variables' files, detect PERMON version
##########################################################
set (permon_conf_rules "${PERMON_DIR}/lib/permon/conf/permon_rules")
set (permon_conf_variables "${PERMON_DIR}/lib/permon/conf/permon_variables")

# All remaining is under this condition except the standard variable handling at the very end
if (permon_conf_rules AND permon_conf_variables AND NOT permon_config_current)
  #permon_get_version()

  # Put variables into environment since they are needed to get
  # configuration (permonvariables) in the PERMON makefile
  set (ENV{PERMON_DIR} "${PERMON_DIR}")
  set (ENV{PERMON_ARCH} "${PETSC_ARCH}")

  ################################################################
  # Extracting  variables from PERMON configuration
  ##############################################################
  # A temporary makefile to probe the PERMON configuration
  set (permon_config_makefile "${PROJECT_BINARY_DIR}/Makefile.permon")
  file (WRITE "${permon_config_makefile}"
"## This file was autogenerated by FindPERMON.cmake
# PERMON_DIR  = ${PERMON_DIR}
# PERMON_ARCH = ${PERMON_ARCH}
include ${permon_conf_rules}
include ${permon_conf_variables}
show :
\t-@echo -n \${\${VARIABLE}}
")

  ######################################
  macro (PERMON_GET_VARIABLE name var)
    set (${var} "NOTFOUND" CACHE INTERNAL "Cleared" FORCE)
    execute_process (COMMAND ${MAKE_EXECUTABLE} --no-print-directory -f ${permon_config_makefile} show VARIABLE=${name}
      OUTPUT_VARIABLE ${var}
      RESULT_VARIABLE permon_return)
  endmacro (PERMON_GET_VARIABLE)

  macro(PERMON_EXPORT_VARIABLES var_list)
      #message(STATUS ${var_list})   
      
      foreach( var ${var_list} )
          permon_get_variable(${var} PERMON_VAR_tmp_${var})
          message(STATUS "EXPORTING PERMON VARIABLE: " ${var} " as PERMON_VAR_${var} = " ${PERMON_VAR_${var}})
          set(PERMON_VAR_${var} ${PERMON_VAR_tmp_${var}} CACHE STRING "Exported variable from PERMON." FORCE)
      endforeach(var)
  endmacro(PERMON_EXPORT_VARIABLES)  
  
  permon_get_variable (PERMON_LIB_DIR           permon_lib_dir)
  permon_get_variable (PERMON_EXTERNAL_LIB      permon_libs_external)
  permon_get_variable (PERMON_CCPPFLAGS         permon_cpp_line)
  permon_get_variable (PERMON_INCLUDE           permon_include)
  permon_get_variable (PERMON_DEP_LIBS          permon_deps)
  permon_get_variable (PCC                      permon_cc)
  permon_get_variable (PCC_FLAGS                permon_cc_flags)
  permon_get_variable (MPIEXEC                  permon_mpiexec)
  
  set(export_variable_list ${PERMON_EXPORT_LIST} ) 
  permon_export_variables("${export_variable_list}") 
  

  # We are done with the temporary Makefile, calling PERMON_GET_VARIABLE after this point is invalid!
  file (REMOVE ${permon_config_makefile})
  # add libraries specified by user (possibly fixing wrong sequence provided by PERMON)
  set(permon_libs_external "${permon_libs_external} ${PERMON_ADDITIONAL_LIBS}")
  

  include (ResolveCompilerPaths)
  # Extract include paths and libraries from compile command line
  resolve_includes (permon_includes_all "${permon_cpp_line}")


  #TODO PERMON was never tested on windows
  #on windows we need to make sure we're linking against the right
  #runtime library
  if (WIN32)
    if (permon_cc_flags MATCHES "-MT")
      set(using_md False)
      foreach(flag_var
          CMAKE_C_FLAGS CMAKE_C_FLAGS_DEBUG CMAKE_C_FLAGS_RELEASE
          CMAKE_C_FLAGS_MINSIZEREL CMAKE_C_FLAGS_RELWITHDEBINFO
          CMAKE_CXX_FLAGS CMAKE_CXX_FLAGS_DEBUG CMAKE_CXX_FLAGS_RELEASE
          CMAKE_CXX_FLAGS_MINSIZEREL CMAKE_CXX_FLAGS_RELWITHDEBINFO)
        if(${flag_var} MATCHES "/MD")
          set(using_md True)
        endif(${flag_var} MATCHES "/MD")
      endforeach(flag_var)
      if(${using_md} MATCHES "True")
        message(WARNING "PERMON was built with /MT, but /MD is currently set.
 See http://www.cmake.org/Wiki/CMake_FAQ#How_can_I_build_my_MSVC_application_with_a_static_runtime.3F")
      endif(${using_md} MATCHES "True")
    endif (permon_cc_flags MATCHES "-MT")
  endif (WIN32)

  include (CorrectWindowsPaths)
  convert_cygwin_path(permon_lib_dir)
  #message (STATUS "permon_lib_dir ${permon_lib_dir}")

  ###############################################################################
  # Find PERMON libraries
  ##############################################################################
  macro (PERMON_FIND_LIBRARY suffix name)
    set (PERMON_LIBRARY_${suffix} "NOTFOUND" CACHE INTERNAL "Cleared" FORCE) # Clear any stale value, if we got here, we need to find it again
    if (WIN32)
      set (libname lib${name}) #windows expects "libfoo", linux expects "foo"
    else (WIN32)
      set (libname ${name})
    endif (WIN32)
    find_library (PERMON_LIBRARY_${suffix} NAMES ${libname} HINTS ${permon_lib_dir} NO_DEFAULT_PATH)
    set (PERMON_LIBRARIES_${suffix} "${PERMON_LIBRARY_${suffix}}")
    mark_as_advanced (PERMON_LIBRARY_${suffix})
  endmacro (PERMON_FIND_LIBRARY suffix name)

  #TODO test with multiple pkgs
  # Look for permonvec first, if it doesn't exist, we must be using single-library
  permon_find_library (SYS permonvec)
  if (PERMON_LIBRARY_SYS)
    permon_find_library (VEC  permonvec)
    permon_find_library (MAT  permonmat)
    permon_find_library (PC   permonpc)
    permon_find_library (KSP  permonksp)
    permon_find_library (QP   permonQP)
    permon_find_library (QPC  permonQPC)
    permon_find_library (QPPF permonQPPF)
    permon_find_library (QPS  permonQPS)
    macro (PERMON_JOIN libs deps)
      list (APPEND PERMON_LIBRARIES_${libs} ${PERMON_LIBRARIES_${deps}})
      #message(STATUS "PERMON_LIBRARIES_${libs}: " ${PERMON_LIBRARIES_${libs}})       
    endmacro (PERMON_JOIN libs deps)
    permon_join (VEC   SYS)
    permon_join (MAT   VEC)
    permon_join (QP  MAT)
    permon_join (QPC QP)
    permon_join (QPPF QPC)
    permon_join (PC QPPF)
    permon_join (QPS PC)
    permon_join (KSP QPS)
    permon_join (ALL  KSP)
  else ()
    set (PERMON_LIBRARY_VEC "NOTFOUND" CACHE INTERNAL "Cleared" FORCE) # There is no libpermonvec
    permon_find_library (SINGLE permon)
    foreach (pkg SYS MAT KSP PC QP QPC QPPF QPS ALL)
      set (PERMON_LIBRARIES_${pkg} "${PERMON_LIBRARY_SINGLE}")
    endforeach ()
  endif ()
  if (PERMON_LIBRARY_QPS)
    message (STATUS "Recognized PERMON install with separate libraries for each package")
  else ()
    message (STATUS "Recognized PERMON install with single library for all packages")
  endif ()

  include(Check${PERMON_LANGUAGE_BINDINGS}SourceRuns)
  include(Check${PERMON_LANGUAGE_BINDINGS}SourceCompiles)
  
  #########################################################################
  # Test that PERMON works
  ############################################################
  macro (PERMON_TEST_RUNS includes libraries runs)
    if(${PERMON_LANGUAGE_BINDINGS} STREQUAL "C")
      set(_PERMON_ERR_FUNC "CHKERRQ(ierr)")
    elseif(${PERMON_LANGUAGE_BINDINGS} STREQUAL "CXX")
      set(_PERMON_ERR_FUNC "CHKERRXX(ierr)")
    endif()

    set(_PERMON_TEST_SOURCE "
static const char help[] = \"PERMON test program.\";
#include <permonqp.h>
int main(int argc,char *argv[]) {
  PetscErrorCode ierr;
  QP qp;

  ierr = PermonInitialize(&argc,&argv,0,help);if (ierr) return ierr;
  ierr = QPCreate(PETSC_COMM_WORLD,&qp);${_PERMON_ERR_FUNC};
  ierr = QPSetFromOptions(qp);${_PERMON_ERR_FUNC};
  ierr = QPDestroy(&qp);${_PERMON_ERR_FUNC};
  ierr = PermonFinalize();${_PERMON_ERR_FUNC};
  return ierr;
}
")

    # check if we can at least compile the source
	  message (STATUS " compiling with ${includes} ${libraries}")

    multipass_c_source_compiles("${includes}" "${libraries}" "${_PERMON_TEST_SOURCE}" ${runs})

    if (${${runs}}) 
      # check if we can run the executable 
      multipass_c_source_runs ("${includes}" "${libraries}" "${_PERMON_TEST_SOURCE}" ${runs}_runs)

      if (${${runs}_runs}) 
        set (PERMON_EXECUTABLE_RUNS "YES" CACHE BOOL
          "Can the system successfully compile and link a PERMON executable?  This variable can be manually set to \"YES\" to force CMake to accept a given PERMON configuration.  If you change PERMON_DIR, or PERMON_CURRENT you would have to reset this variable." FORCE)         
      endif(${${runs}_runs})    
    endif (${${runs}}) 
  endmacro (PERMON_TEST_RUNS)


  find_path (PERMON_INCLUDE_DIR permonqp.h HINTS "${PERMON_DIR}" PATH_SUFFIXES include NO_DEFAULT_PATH)
  mark_as_advanced (PERMON_INCLUDE_DIR)
  set (permon_includes_minimal ${PERMON_INCLUDE_CONF} ${PERMON_INCLUDE_DIR})

  # Macro resolve_libraries comes from ResolveCompilerPaths, it tries resolve all libraries form given compiler line
  # Unfortunately this only mimics compiler resolutions and occasionally can be incorrect. In such a case use PERMON_ADDITIONAL_LIBS.
  #message(STATUS "[FindPERMON] Try to resolve libraries from: '${permon_libs_external}'")
  resolve_libraries (permon_libraries_external "${permon_libs_external}")
  #message(STATUS "[FindPERMON] Resolved path: '${permon_libraries_external}'")

  #message(STATUS "[FindPERMON] Try to resolve libraries from: '${permon_deps}'")
  resolve_libraries (permon_libraries_deps "${permon_deps}")
  #message(STATUS "[FindPERMON] Resolved path: '${permon_libraries_deps}'")

  set (permon_libs_needed ${PERMON_LIBRARIES_QPS})
  
  foreach (pkg SYS MAT QP QPC QPPF PC QPS KSP ALL)
    list (APPEND PERMON_LIBRARIES_${pkg}  ${permon_libraries_external})
    list (APPEND PERMON_LIBRARIES_${pkg}  ${permon_libraries_deps})
  endforeach (pkg)

  message(STATUS "lib qps: '${PERMON_LIBRARIES_QPS}'")
  permon_test_runs ("${permon_includes_minimal}" "${PERMON_LIBRARIES_QPS}" permon_works_alllibraries)
  if (permon_works_alllibraries)
     message (STATUS "PERMON only need minimal includes, but requires explicit linking to all dependencies.  This is expected when PERMON is built with static libraries.")
    set (permon_includes_needed ${permon_includes_minimal})
  else (permon_works_alllibraries)
  
    # Multipass_test_4 ####################      
    # It looks like we really need everything, should have listened to Matt
    set (permon_includes_needed ${permon_includes_all})
    permon_test_runs ("${permon_includes_all}" "${PERMON_LIBRARIES_QPS}" permon_works_all)
    if (permon_works_all) # We fail anyways
      message (STATUS "PERMON requires extra include paths and explicit linking to all dependencies.  This probably means you have static libraries and something unexpected in PERMON headers.")
	    set (permon_includes_needed ${permon_includes_all})
    else (permon_works_all) # We fail anyways
	    message (STATUS "
	Can not compile and link PERMON executable, probably some missing libraries.
	See BUILD_DIR/CMakeFiles/CMakeError.log for reasons, check library resolution after 'MULTIPASS_TEST_*_permon_works_allincludes'.
	Try to add missing libraries through PERMON_ADDITIONAL_LIBS variable.")
    endif (permon_works_all)
  endif (permon_works_alllibraries)
  
  if (permon_includes_needed)            # this indicates PERMON_FOUND
    if (${PERMON_EXECUTABLE_RUNS})       # this is optional
    else()
        message(STATUS "
        Can compile but can not run PERMON executable. Probably problem with mpiexec or with shared libraries.
        See See BUILD_DIR/CMakeFiles/CMakeError.log for reasons.")
    endif(${PERMON_EXECUTABLE_RUNS})    
  endif(permon_includes_needed)  

  # reset to minimal
  # We do an out-of-source build so __FILE__ will be an absolute path, hence __INSDIR__ is superfluous
  set (PERMON_DEFINITIONS "-D__INSDIR__=" CACHE STRING "PERMON definitions" FORCE)
  # Sometimes this can be used to assist FindMPI.cmake
  #set (PERMON_MPIEXEC ${permon_mpiexec} CACHE FILEPATH "Executable for running PERMON MPI programs" FORCE)
  set (PERMON_INCLUDES ${permon_includes_minimal} CACHE STRING "PERMON include path" FORCE)
  set (PERMON_LIBRARY ${permon_libs_needed} CACHE STRING "PERMON libraries" FORCE)
  #set (PERMON_COMPILER ${permon_cc} CACHE FILEPATH "PERMON compiler" FORCE)
  # TODO do not ignore potential external libs
  #set (PERMON_EXTERNAL_LIB ${permon_libraries_deps} CACHE STRING "PERMON external libraries" FORCE)
  set (PERMON_EXTERNAL_LIB ${permon_libraries_external} CACHE STRING "PERMON external libraries" FORCE)
  # Note that we have forced values for all these choices.  If you
  # change these, you are telling the system to trust you that they
  # work.  It is likely that you will end up with a broken build.
  mark_as_advanced (PERMON_INCLUDES PERMON_LIBRARIES PERMON_DEFINITIONS PERMON_EXECUTABLE_RUNS)
endif ()

include (FindPackageHandleStandardArgs)
find_package_handle_standard_args (PERMON
  "PERMON could not be found.  Be sure to set PERMON_DIR"
  PERMON_INCLUDES PERMON_LIBRARY)
