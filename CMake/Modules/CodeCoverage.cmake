# Edited https://github.com/bilke/cmake-modules/blob/master/CodeCoverage.cmake
# This library can be used to collect code coverage and generate HTML or XML report
# This library uses standard CMAKE variables:
# - determining CMAKE compiler  (CMAKE_COMPILER_IS_GNUCXX and CMAKE_CXX_COMPILER_ID)
# - resolving source dir        (CMAKE_SOURCE_DIR)
# - resolving binary dir        (CMAKE_BINARY_DIR)
#
# note: no CMAKE variables are affected
#
# usage:
# 1) call function in your cmakelists which adds custom target (don't forget to include this file)
#    specify custom target name (_targetname) and XML/HTML report name (_outputname)
#    example: 
#       INCLUDE (CodeCoverage)
#       COVERAGE_SETUP_COBERTURA_COLLECTION (collect-coverage coverage.xml)
# 2) compile source code with flags supporting code coverage AND without any
#    optimization (optimization can alter source code in some cases)
#    example flags: -g -O0 -fprofile-arcs -ftest-coverage
# 3) run your program, run tests
# 4) collect coverage
#    example: 
#       make collect-coverage
#
# coverage.xml should contain cobertura-style XML coverage


# Check prereqs
FIND_PROGRAM( GCOV_PATH gcov )
FIND_PROGRAM( LCOV_PATH lcov )
FIND_PROGRAM( GENHTML_PATH genhtml )
FIND_PROGRAM( GCOVR_PATH gcovr PATHS ${CMAKE_SOURCE_DIR}/tests)

IF(NOT GCOV_PATH)
	MESSAGE(STATUS "gcov not found! Aborting...")
ENDIF() # NOT GCOV_PATH

IF(NOT CMAKE_COMPILER_IS_GNUCXX)
	# Clang version 3.0.0 and greater now supports gcov as well.
	MESSAGE(WARNING "Compiler is not GNU gcc! Clang Version 3.0.0 and greater supports gcov as well, but older versions don't.")
	
	IF(NOT "${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang")
		MESSAGE(FATAL_ERROR "Compiler is not GNU gcc! Aborting...")
	ENDIF()
ENDIF() # NOT CMAKE_COMPILER_IS_GNUCXX




# Function creates HTML report by processing code coverage files
#
# Param _targetname     The name of new the custom make target
# Param _outputname     lcov output is generated as _outputname.info
#                       HTML report is generated in _outputname/index.html
FUNCTION(COVERAGE_SETUP_HTML_COLLECTION _targetname _outputname)

	IF(NOT LCOV_PATH)
		MESSAGE(FATAL_ERROR "lcov not found! Cannot generate HTML coverage")
	ENDIF() # NOT LCOV_PATH

	IF(NOT GENHTML_PATH)
		MESSAGE(FATAL_ERROR "genhtml not found! Cannot generate HTML coverage")
	ENDIF() # NOT GENHTML_PATH

	# Setup target
	ADD_CUSTOM_TARGET(${_targetname}
		
		# Cleanup lcov
		${LCOV_PATH} --directory . --zerocounters
		
		# Capturing lcov counters and generating report
		COMMAND ${LCOV_PATH} --directory . --capture --output-file ${_outputname}.info
		COMMAND ${LCOV_PATH} --remove ${_outputname}.info 'tests/*' '/usr/*' --output-file ${_outputname}.info.cleaned
		COMMAND ${GENHTML_PATH} -o ${_outputname} ${_outputname}.info.cleaned
		COMMAND ${CMAKE_COMMAND} -E remove ${_outputname}.info ${_outputname}.info.cleaned
		
		WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
		COMMENT "Resetting code coverage counters to zero.\nProcessing code coverage counters and generating report."
	)
	
	# Show info where to find the report
	ADD_CUSTOM_COMMAND(TARGET ${_targetname} POST_BUILD
		COMMAND ;
		COMMENT "Open ./${_outputname}/index.html in your browser to view the coverage report."
	)

ENDFUNCTION() # COVERAGE_SETUP_HTML_COLLECTION




# Function creates XML report by processing code coverage files
#
# Param _targetname     The name of new the custom make target
# Param _outputname     cobertura output is generated as _outputname (e.g. coverage.xml)
FUNCTION(COVERAGE_SETUP_COBERTURA_COLLECTION _targetname _outputname)

	IF(NOT GCOVR_PATH)
		MESSAGE(FATAL_ERROR "gcovr not found! Aborting...")
	ENDIF() # NOT GCOVR_PATH

	ADD_CUSTOM_TARGET(${_targetname}

		# Running gcovr
		COMMAND ${GCOVR_PATH} -x --xml-pretty -r ${CMAKE_SOURCE_DIR} -k -v -o ${_outputname} -s -e "\\.\\*\\\\\\.hpp$$" -e "\\.\\*\\\\\\.h$$" -e "\\.\\*\\\\\\.hh$$"
		WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
		COMMENT "Running gcovr to produce Cobertura code coverage report."
	)

	# Show info where to find the report
	ADD_CUSTOM_COMMAND(TARGET ${_targetname} POST_BUILD
		COMMAND ;
		COMMENT "Cobertura code coverage report saved in ${_outputname}."
	)

ENDFUNCTION() # COVERAGE_SETUP_COBERTURA_COLLECTION
