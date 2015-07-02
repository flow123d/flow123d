# Locate yaml-cpp
#
# This module defines
#  YAMLCPP_FOUND, if false, do not try to link to yaml-cpp
#  YAMLCPP_LIBRARY, where to find yaml-cpp
#  YAMLCPP_INCLUDE_DIR, where to find yaml.h
#
# By default, the dynamic libraries of yaml-cpp will be found. To find the static ones instead,
# you must set the YamlCpp_STATIC_LIBRARY variable to TRUE before calling find_package(YamlCpp ...).
#
# If yaml-cpp is not installed in a standard path, you can use the YamlCpp_DIR CMake variable
# to tell CMake where yaml-cpp is.

# attempt to find static library first if this is set
# currently force static library
#if(YamlCpp_STATIC_LIBRARY)
    set(YamlCpp_STATIC libyaml-cpp.a)
#endif()

# find the yaml-cpp include directory
find_path(YamlCpp_INCLUDE_DIR yaml-cpp/yaml.h
          PATH_SUFFIXES include
          PATHS
          /usr/local/include/
          /usr/include/
          ${YamlCpp_DIR}/include/)

# find the yaml-cpp library
find_library(YamlCpp_LIBRARY
             NAMES ${YamlCpp_STATIC} yaml-cpp
             PATH_SUFFIXES lib64 lib
             PATHS  /usr/local
                    /usr
                    ${YamlCpp_DIR}/lib)

# handle the QUIETLY and REQUIRED arguments and set YAMLCPP_FOUND to TRUE if all listed variables are TRUE
include(FindPackageHandleStandardArgs)
FIND_PACKAGE_HANDLE_STANDARD_ARGS(YamlCpp DEFAULT_MSG YamlCpp_INCLUDE_DIR YamlCpp_LIBRARY)
mark_as_advanced(YamlCpp_INCLUDE_DIR YamlCpp_LIBRARY)
