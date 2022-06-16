# Find PERMON library
# Input: PERMON_ROOT
#
# Output:
# PERMON_INCLUDES - where to find header files
# PERMON_LIBRARIES - The link sequence for library and dependencies.
# PERMON_FOUND - True if found.


########################################################## 
# Use package multipass to alow multiple tries when search for library.
#
# FIND_PACKAGE_MULTIPASS (Name CURRENT
#  STATES VAR0 VAR1 ...
#  DEPENDENTS DEP0 DEP1 ...)
#
#  This function creates a cache entry <UPPERCASED-Name>_CURRENT which
#  the user can set to "NO" to trigger a reconfiguration of the package.
#  The first time this function is called, the values of
#  <UPPERCASED-Name>_VAR0, ... are saved.  If <UPPERCASED-Name>_CURRENT
#  is false or if any STATE has changed since the last time
#  FIND_PACKAGE_MULTIPASS() was called, then CURRENT will be set to "NO",
#  otherwise CURRENT will be "YES".  IF not CURRENT, then
#  <UPPERCASED-Name>_DEP0, ... will be FORCED to NOTFOUND.
##########################################################
include (FindPackageMultipass)
find_package_multipass (PERMON config_current
  STATES ROOT
  DEPENDENTS LAST_ROOT INCLUDES LIBRARIES FOUND)


message(STATUS "PERMON_ROOT: ${PERMON_ROOT}")

find_path (PERMON_INCLUDES permonksp.h  HINTS ${PERMON_ROOT}/include )
find_library (PERMON_LIBRARY NAMES permon  HINTS ${PERMON_ROOT}/lib )

if (PERMON_INCLUDES AND PERMON_LIBRARY)
  message(STATUS "found, include: ${PERMON_INCLUDES}, lib: ${PERMON_LIBRARY} ")
  set(PERMON_LAST_ROOT ${PERMON_ROOT} CACHE FILEPATH "_Clearad_") 
  #string(REPLACE REGEX "/src$" "" PERMON_ROOT ${PERMON_INCLUDES})
  #find_program (MAKE_EXECUTABLE NAMES make gmake)

  set(PERMON_LIBRARIES ${PERMON_LIBRARY})

#     ###################################
#     macro (PERMON_GET_VARIABLE name var)
#       set (${var} "${var}_NOTFOUND" CACHE INTERNAL "Cleared" FORCE)
#       execute_process (COMMAND ${MAKE_EXECUTABLE} --no-print-directory -f ${config_makefile} show VARIABLE=${name}
#         OUTPUT_VARIABLE ${var}
#         RESULT_VARIABLE make_return)
#       message(STATUS " exporting ${name} into cmake variable ${var} =  " ${${var}} )
#     endmacro (PERMON_GET_VARIABLE)
# 
#     bddcml_get_variable("PERMON_LINK_SEQUENCE" PERMON_LINK_SEQ)
#     bddcml_get_variable("PERMON_MPCFLAGS" PERMON_CDEFS)
#     bddcml_get_variable("PERMON_INC" PERMON_INC)
#     STRING(REGEX REPLACE "( |^)-I" ";"  PERMON_INC_DIRS ${PERMON_INC})
#     message(STATUS ${PERMON_INC_DIRS})
#     STRING(REGEX REPLACE " *" ";"  PERMON_INC_LIST ${PERMON_INC_DIRS}) 
#     set(PERMON_INCLUDES ${PERMON_INCLUDES} ${PERMON_INC_DIRS})
#     
#   
#     message(STATUS "PERMON_LIBS: ${PERMON_LINK_SEQ}")
#     message(STATUS "PERMON_CDEFS: ${PERMON_CDEFS}")
# 
#     file (REMOVE ${config_makefile})
#   elseif()
#     set(PERMON_INCLUDES "")
#     set(PERMON_LIBRARY "")
#   endif()  
endif()  



#include (ResolveCompilerPaths)
# Extract include paths and libraries from compile command line
#resolve_LIBRARIES(PERMON_LIBRARIES "${PERMON_LINK_SEQ}")
#message(STATUS "B SEQ:" ${PERMON_LINK_SEQ})
# Extract BLOPEX object files 
#STRING(REGEX MATCHALL "[^ ]*\\.o" BLOPEX_OBJS "${PERMON_LINK_SEQ}")
#message(STATUS "BOBJ: ${BLOPEX_OBJS}")
#STRING(REGEX REPLACE <regular_expression>  <replace_expression> <output variable>  <input> [<input>...])
# ATTENTION: input has to be quoted otherwise it interprets contens of variable, in particular run REGEX on every item in ; separated list
#STRING(REGEX REPLACE "([^;]*BLOPEX\\.a)" "${BLOPEX_OBJS};\\1" PERMON_LIBRARIES "${PERMON_LIBS}" )
#STRING(REGEX MATCH "[^;]*BLOPEX\\.a"  PERMON_BLOPEX_LIB "${PERMON_LIBS}" )
#STRING(REGEX MATCH "[^;]"  PERMON_BEGIN_LIB ${PERMON_LIBS} )
#message(STATUS "out: ${PERMON_LIBRARIES}")
#message(STATUS "out: ${PERMON_BLOPEX_LIB}")
#message(STATUS "in: ${PERMON_LIBS}")


message(STATUS "found: ${PERMON_INCLUDES}, ${PERMON_LIBRARIES}")


# handle the QUIETLY and REQUIRED arguments and set FFTW_FOUND to TRUE if
# all listed variables are TRUE
include (FindPackageHandleStandardArgs)
find_package_handle_standard_args (PERMON DEFAULT_MSG PERMON_LIBRARIES PERMON_INCLUDES)

mark_as_advanced (PERMON_LIBRARIES PERMON_INCLUDES)



