# Locate zlib
#
# This module defines
#  ZLIB_FOUND, if false, do not try to link to zlib
#  ZLIB_LIBRARY, where to find zlib
#  ZLIB_INCLUDE_DIR, where to find ZLibrary.h
#
# By default, the dynamic libraries of zlib will be found. To find the static ones instead,
# you must set the Zlib_STATIC_LIBRARY variable to TRUE before calling find_package(Zlib ...).
#
# If zlib is not installed in a standard path, you can use the Zlib_ROOT_HINT CMake variable
# to tell CMake where zlib is.

# attempt to find static library first if this is set
# currently force static library

#if(Zlib_STATIC_LIBRARY)
    set(Zlib_STATIC libz.so.1)
#endif()

# find the zlib include directory
find_path(Zlib_INCLUDE_DIR zlib.h
#          PATH_SUFFIXES zlibrary
#          HINTS 
#              ${Zlib_ROOT_HINT}/include/
#          PATHS
#              /usr/local/include/
#              /usr/include/
          )

# find the zlib library
find_library(Zlib_LIBRARY
             NAMES ${Zlib_STATIC} Zlib
#             PATH_SUFFIXES lib64 lib
#             HINTS
#                 ${Zlib_ROOT_HINT}/lib/
#             PATHS  
#                /usr/local
#                /usr
             )

# handle the QUIETLY and REQUIRED arguments and set ZLIB_FOUND to TRUE if all listed variables are TRUE
include(FindPackageHandleStandardArgs)
FIND_PACKAGE_HANDLE_STANDARD_ARGS(Zlib DEFAULT_MSG Zlib_INCLUDE_DIR Zlib_LIBRARY)
mark_as_advanced(Zlib_INCLUDE_DIR Zlib_LIBRARY)
