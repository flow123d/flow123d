# Locate Gmsh 4.0
#
# This module defines
#  GMSH_FOUND, if false, do not try to link to gmsh
#  GMSH_LIBRARY, where to find gmsh
#  GMSH_INCLUDE_DIR, where to find gmsh.h
#
# By default, the dynamic libraries of gmsh will be found. To find the static ones instead,
# you must set the Gmsh_STATIC_LIBRARY variable to TRUE before calling find_package(Gmsh ...).
#
# If Gmsh is not installed in a standard path, you can use the Gmsh_ROOT_HINT CMake variable
# to tell CMake where gmsh is.

# attempt to find static library first if this is set
# currently force static library

#if(Gmsh_STATIC_LIBRARY)
    set(Gmsh_STATIC libgmsh.so.4.9)
#endif()

# find the Gmsh include directory
find_path(Gmsh_INCLUDE_DIR gmsh.h
          PATH_SUFFIXES include
          HINTS 
              ${Gmsh_ROOT_HINT}/include/
          PATHS
              /usr/local/include/
              /usr/include/
          )

# find the Gmsh library
find_library(Gmsh_LIBRARY
             NAMES ${Gmsh_STATIC} gmsh
             PATH_SUFFIXES lib64 lib
             HINTS
                 ${Gmah_ROOT_HINT}/lib/
             PATHS  
                /usr/local
                /usr
             )

# handle the QUIETLY and REQUIRED arguments and set GMSH_FOUND to TRUE if all listed variables are TRUE
include(FindPackageHandleStandardArgs)
FIND_PACKAGE_HANDLE_STANDARD_ARGS(Gmsh DEFAULT_MSG Gmsh_INCLUDE_DIR Gmsh_LIBRARY)
mark_as_advanced(Gmsh_INCLUDE_DIR Gmsh_LIBRARY)
