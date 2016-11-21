add_test_directory( "" )

# in order to use these tests additional files will be copied to CMAKE_CURRENT_BINARY_DIR
# we also add target for compiling consumer.cc tester-like c++ program which simulates 
# memory and duration stress
execute_process(COMMAND
    ${CMAKE_COMMAND} -E copy_directory ${CMAKE_CURRENT_SOURCE_DIR}/extras/ ${CMAKE_CURRENT_BINARY_DIR}/extras/)
add_executable(consumer_cc EXCLUDE_FROM_ALL extras/consumer.cc) 


# find python 2.7
find_package(PythonInterp 2.7 REQUIRED)

# Custom target which runs all tests for python scripts. By default output is hidden (flag -b)
# but on error it is displayed (sufficient for all applications)
add_custom_target(test_scripts
    COMMAND ${PYTHON_EXECUTABLE} -m unittest discover -s ${CMAKE_CURRENT_SOURCE_DIR} -bv
    DEPENDS consumer_cc
    WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})

# Custom target which runs all tests for python scripts where output is shown. Since unittests
# are testing correct and also incorrect behaviour output can be misleading (exceptions are thrown
# and displayed yet test succeeds)
add_custom_target(test_scripts_debug
    COMMAND ${PYTHON_EXECUTABLE} -m unittest discover -s ${CMAKE_CURRENT_SOURCE_DIR} -v
WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})