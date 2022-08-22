# Locate pugi-xml
#
# This module defines
#  PUGIXML_FOUND, if false, do not try to link to pugixml
#  PUGIXML_LIBRARY, where to find pugixml
#  PUGIXML_INCLUDE_DIR, where to find pugixml.hpp
#
# By default, the dynamic libraries of pugi-xml will be found. To find the static ones instead,
# you must set the PugiXml_STATIC_LIBRARY variable to TRUE before calling find_package(PugiXml ...).
#
# If pugi-xml is not installed in a standard path, you can use the PugiXml_ROOT_HINT CMake variable
# to tell CMake where pugi-xml is.

# attempt to find static library first if this is set
# currently force static library

#if(PugiXml_STATIC_LIBRARY)
    set(PugiXml_STATIC libpugixml.so.1.10)
#endif()

# find the pugi-xml include directory
find_path(PugiXml_INCLUDE_DIR pugixml.hpp
          PATH_SUFFIXES include
          HINTS 
              ${PugiXml_ROOT_HINT}/include/
          PATHS
              /usr/local/include/
              /usr/include/
          )

# find the pugi-xml library
find_library(PugiXml_LIBRARY
             NAMES ${PugiXml_STATIC} pugixml
             PATH_SUFFIXES lib64 lib
             HINTS
                 ${PugiXml_ROOT_HINT}/lib/
             PATHS  
                /usr/local
                /usr
             )

# handle the QUIETLY and REQUIRED arguments and set PUGIXML_FOUND to TRUE if all listed variables are TRUE
include(FindPackageHandleStandardArgs)
FIND_PACKAGE_HANDLE_STANDARD_ARGS(PugiXml DEFAULT_MSG PugiXml_INCLUDE_DIR PugiXml_LIBRARY)
mark_as_advanced(PugiXml_INCLUDE_DIR PugiXml_LIBRARY)
