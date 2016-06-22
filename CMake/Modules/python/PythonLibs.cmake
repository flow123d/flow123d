# PythonLibs module
#  - this module find PythonLibs package and sets up variables:
#    - PYTHONLIBS_VERSION_MAJOR
#           - major version in python 2.7.1 it'll be number 2
#    - PYTHONLIBS_VERSION_MINOR
#           - minor version in python 2.7.1 it'll be number 7
#    - PYTHONLIBS_VERSION_PATCH
#           - pacth version in python 2.7.1 it'll be number 1

# find package or die
find_package(PythonLibs 2.7 REQUIRED)

include_directories(${PYTHON_INCLUDE_DIRS})
flow_define(HAVE_PYTHON)
# parse version
string(REGEX REPLACE "([0-9]+).*" "\\1" PYTHONLIBS_VERSION_MAJOR ${PYTHONLIBS_VERSION_STRING})
string(REGEX REPLACE "[0-9]+.([0-9]+).*" "\\1" PYTHONLIBS_VERSION_MINOR ${PYTHONLIBS_VERSION_STRING})
string(REGEX REPLACE "[0-9]+.[0-9]+.([0-9]+).*" "\\1" PYTHONLIBS_VERSION_PATCH ${PYTHONLIBS_VERSION_STRING})

# add defines
flow_define_constant(PYTHONLIBS_VERSION_MAJOR ${PYTHONLIBS_VERSION_MAJOR})
flow_define_constant(PYTHONLIBS_VERSION_MINOR ${PYTHONLIBS_VERSION_MINOR})
flow_define_constant(PYTHONLIBS_VERSION_PATCH ${PYTHONLIBS_VERSION_PATCH})
flow_define(PYTHONLIBS_VERSION_STRING ${PYTHONLIBS_VERSION_STRING})

# debug info
message(STATUS "Python version:       ${PYTHONLIBS_VERSION_STRING}")
message(STATUS "Python major version: ${PYTHONLIBS_VERSION_MAJOR}")
message(STATUS "Python minor version: ${PYTHONLIBS_VERSION_MINOR}")
message(STATUS "Python patch version: ${PYTHONLIBS_VERSION_PATCH}")
message(STATUS "Python include:       ${PYTHON_INCLUDE_DIRS}")
message(STATUS "Python libs:          ${PYTHON_LIBRARIES}")