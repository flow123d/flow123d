# MACRO for testing presence of a timer
#   timer_name              name of the variable which will holds result
#   timer_script_location   path to the cpp file which test presence of a timer
#                           and prints its resolution in microseconds
MACRO(TEST_TIMER timer_name timer_script_location)
    SET(TIMER_RESOLUTION No)
    try_run(
        TIMER_EXIT "${timer_name}"
        ${CMAKE_BINARY_DIR} "${timer_script_location}"
        CMAKE_FLAGS "${CMAKE_CXX_FLAGS}"
        RUN_OUTPUT_VARIABLE TIMER_RESOLUTION COMPILE_OUTPUT_VARIABLE TIMER_COMPILE) 
    MESSAGE(STATUS "${timer_name}:\n ${${timer_name}} ${TIMER_RESOLUTION} [us]\n\n")
    if (NOT ${timer_name})
        MESSAGE (STATUS "Compilation error: \n${TIMER_COMPILE}")
    endif()
ENDMACRO(TEST_TIMER timer_name timer_script_location)