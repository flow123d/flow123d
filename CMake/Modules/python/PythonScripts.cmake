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

# setup some variables
set(PYTHON_3RD_PARTY ${PY_BUILD_PREFIX}/lib/python2.7/site-packages)
set(PYTHON_RUNTEST   ${CMAKE_SOURCE_DIR}/src/py123d/bin/runtest.py)

# macro for installing python modules using pip
#    package_name - name of the module to be installed
#    package_dir  - name of the folder which should be checked whether exists 
#                   some packages name and dir are different e.g. pyyaml -> yaml
macro(install_python_lib package_name package_dir)
    if(EXISTS ${PYTHON_3RD_PARTY}/${package_dir})
        message(STATUS "Installing ${package_name} with pip .... skipped, already exists")
    else()
        message(STATUS "Installing ${package_name} with pip ....")
        message(STATUS "-- cd ${Flow123d_SOURCE_DIR}/CMake/Modules/python && pip install --target=${PYTHON_3RD_PARTY} ${package_name}")
        # we start pip install from directory where setup.cfg is placed and overwrite 
        # default python pip configuration, thus avoiding issue with --target usage
        # https://github.com/pypa/pip/pull/3450 and https://github.com/pypa/pip/issues/3056
        # for now we try to install package using this fix (platform centos)
        execute_process(
            COMMAND pip install --target=${PYTHON_3RD_PARTY} ${package_name}
            WORKING_DIRECTORY ${Flow123d_SOURCE_DIR}/CMake/Modules/python
            RESULT_VARIABLE RETURNCODE
            OUTPUT_FILE ${CMAKE_BINARY_DIR}/python_install.log
            ERROR_FILE ${CMAKE_BINARY_DIR}/python_install.log)
        
        # if installation fails we run pip from current dir and do not override 
        # any configuration
        if(NOT RETURNCODE EQUAL 0)
            execute_process(
                COMMAND pip install --target=${PYTHON_3RD_PARTY} ${package_name}
                RESULT_VARIABLE RETURNCODE
                OUTPUT_FILE ${CMAKE_BINARY_DIR}/python_install.log
                ERROR_FILE ${CMAKE_BINARY_DIR}/python_install.log)
        endif()
        # could not install package
        if(NOT RETURNCODE EQUAL 0)
            message(WARNING "Could not install python package ${package_name}")
        endif()
    endif()
endmacro(install_python_lib)

# install python dependencies
install_python_lib(markdown markdown)
install_python_lib(pyyaml yaml)
install_python_lib(psutil psutil)

# create Python runtest wrapper in src/bin
# also create symlink of python runtest wrapper in tests
message(STATUS "Creating python runtest wrapper")
set(PYTHON_RUNTEST_WRAPPER ${CMAKE_SOURCE_DIR}/bin/runtest.sh)
configure_file(${CMAKE_SOURCE_DIR}/CMake/unix_runtest_template ${PYTHON_RUNTEST_WRAPPER} @ONLY)
if (COPY_INSTEAD_OF_SYMLINK)
	execute_process(COMMAND
	    ${CMAKE_COMMAND} -E copy ${PYTHON_RUNTEST_WRAPPER} ${CMAKE_SOURCE_DIR}/tests/runtest.sh
	)
	execute_process(COMMAND
	    ${CMAKE_COMMAND} -E copy_directory ${PYTHON_RUNTEST_WRAPPER} ${CMAKE_SOURCE_DIR}/tests/runtest.sh
	)
else()
	execute_process(COMMAND
	    ${CMAKE_COMMAND} -E create_symlink ${PYTHON_RUNTEST_WRAPPER} ${CMAKE_SOURCE_DIR}/tests/runtest.sh
	)
endif()

# construct PY_WRAPPER_PATHS for wrapper pythonenv
set(PY_WRAPPER_PATHS "/lib/python2.7")
list(APPEND PY_WRAPPER_PATHS ${PYTHON_EXTRA_MODULES_PATH})
foreach(PY_PATH ${PYTHON_SYSPATH})
    list(APPEND PY_WRAPPER_PATHS "/lib/python2.7/${PY_PATH}")
endforeach()

# configure pythonenv.sh script (populates py wrapper path used in wrapper file)
message(STATUS "Creating pythonenv.sh wrapper")
configure_file(${CMAKE_SOURCE_DIR}/CMake/pythonenv_template ${CMAKE_SOURCE_DIR}/bin/pythonenv.sh @ONLY)