# PythonCopy module
#  - this module will install scripts dependencies and
#    copy python scripts to appropriate folders
#    
#    
#  - Following variables MUST BE SET:
#    - PY_BUILD_PREFIX
#           - directory where libs will be installed
#             in ${PY_BUILD_PREFIX}/lib/python2.7/site-packages


# check working python 2.7+
find_package(PythonInterp 2.7 REQUIRED)

# macro for installing python modules using pip
#    package_name - name of the module to be installed
#    package_dir  - name of the folder which should be checked whether exists 
#                   some packages name and dir are different e.g. pyyaml -> yaml
macro(install_python_lib package_name package_dir)
    execute_process(COMMAND
        pwd
    )
    if(EXISTS ${PYTHON_3RD_PARTY}/${package_dir})
        message(STATUS "Installing ${package_name} with pip .... skipped, already exists")
    else()
        message(STATUS "Installing ${package_name} with pip ....")
        message(STATUS "-- cd ${Flow123d_SOURCE_DIR}/CMake/Modules/python && pip install --target=${PYTHON_3RD_PARTY} ${package_name}")
        execute_process(COMMAND cd ${Flow123d_SOURCE_DIR}/CMake/Modules/python && pip install --target=${PYTHON_3RD_PARTY} ${package_name})
    endif()
endmacro(install_python_lib)

# setup some variables
set(PYTHON_3RD_PARTY ${PY_BUILD_PREFIX}/lib/python2.7/site-packages)
set(PYTHON_RUNTEST   ${CMAKE_SOURCE_DIR}/bin/python/runtest.py)

# install python dependencies
install_python_lib(pyyaml yaml)
install_python_lib(psutil psutil)
install_python_lib(markdown markdown)

# create Python runtest wrapper in src/bin
# also create symlink of python runtest wrapper in tests
message(STATUS "Creating python runtest wrapper")
set(PYTHON_RUNTEST_WRAPPER ${CMAKE_SOURCE_DIR}/bin/runtest.sh)
configure_file(${CMAKE_SOURCE_DIR}/CMake/unix_runtest_template ${PYTHON_RUNTEST_WRAPPER} @ONLY)
execute_process(COMMAND
    ${CMAKE_COMMAND} -E create_symlink ${PYTHON_RUNTEST_WRAPPER} ${CMAKE_SOURCE_DIR}/tests/runtest.sh
)