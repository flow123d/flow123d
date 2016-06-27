# PythonCopy module
#  - this module will install scripts dependencies and
#    copy python scripts to appropriate folders
#    
#    
#  - Following variables MUST BE SET:
#    - PY_BUILD_PREFIX
#           - list of all path which should be added to PYTHONPATH after python
#             is copied out relative to CMAKE_BINARY_DIR
#    - PY_BUILD_PREFIX
#           - directory where 3rd party will be installed
#             in ${PY_BUILD_PREFIX}/lib/python2.7/site-packages


# check working python 2.7+
find_package(PythonInterp 2.7 REQUIRED)

# macro for installing python modules using pip
#    package_name - name of the module to be installed
#    package_dir  - name of the folder which should be checked whether exists 
#                   some packages name and dir are different e.g. pyyaml -> yaml
macro(install_python_lib package_name package_dir)
    if(EXISTS ${PYTHON_3RD_PARTY}/${package_dir})
        message(STATUS "Installing ${package_name} with pip .... skipped, already exists")
    else()
        message(STATUS "Installing ${package_name} with pip ....")
        execute_process(COMMAND pip install --target=${PYTHON_3RD_PARTY} ${package_name})
    endif()
endmacro(install_python_lib)

# setup some variables
set(PYTHON_3RD_PARTY ${PY_BUILD_PREFIX}/lib/python2.7/site-packages)
set(PYTHON_RUNTEST   ${CMAKE_SOURCE_DIR}/bin/python/runtest.py)


# install python dependencies
install_python_lib(pyyaml yaml)
install_python_lib(psutil psutil)
install_python_lib(markdown markdown)


# get current PYTHONPATH variable
execute_process(COMMAND
    ${PYTHON_EXECUTABLE} ${CMAKE_SOURCE_DIR}/bin/python/python_path_script.py sys_path --format cmake
    OUTPUT_VARIABLE PYTHON_PATHS
    OUTPUT_STRIP_TRAILING_WHITESPACE)

# prepare path for wrappers (variable PY_WRAPPER_PATHS)
# everytime we use lib root and 3rd party directory
set(PY_WRAPPER_PATHS
    "/lib/python2.7"
    "/lib/python2.7/site-packages"
    ${PYTHON_EXTRA_MODULES_PATH})

# go through all paths in python sys.path and those in PYTHON_ROOT "transfer" 
# to our lib dir
foreach(ITEM ${PYTHON_PATHS})
    string(FIND ${ITEM} ${PYTHON_ROOT}/ _index)
    if(${_index} EQUAL 0)
        string(REPLACE "${PYTHON_ROOT}/" "" ITEM_NEWPATH ${ITEM})
        list(APPEND PY_WRAPPER_PATHS "/lib/python2.7/${ITEM_NEWPATH}")
    endif()
endforeach()


# create Python runtest wrapper in src/bin
# also create symlink of python runtest wrapper in tests
message(STATUS "Creating python runtest wrapper")
set(PYTHON_RUNTEST_WRAPPER ${CMAKE_SOURCE_DIR}/bin/runtest.sh)
configure_file(${CMAKE_SOURCE_DIR}/CMake/unix_runtest_template ${PYTHON_RUNTEST_WRAPPER} @ONLY)
execute_process(COMMAND
    ${CMAKE_COMMAND} -E create_symlink ${PYTHON_RUNTEST_WRAPPER} ${CMAKE_SOURCE_DIR}/tests/runtest.sh
)

# configure pythonenv.sh script (populates py wrapper path used in wrapper file)
configure_file(${CMAKE_SOURCE_DIR}/CMake/pythonenv_template ${CMAKE_SOURCE_DIR}/bin/pythonenv.sh @ONLY)