# Edited https://github.com/bilke/cmake-modules/blob/master/CodeCoverage.cmake

# Check prereqs
FIND_PROGRAM( GCOV_PATH gcov )
FIND_PROGRAM( LCOV_PATH lcov )
FIND_PROGRAM( GENHTML_PATH genhtml )
FIND_PROGRAM( GCOVR_PATH gcovr PATHS ${CMAKE_SOURCE_DIR}/tests)

IF(NOT GCOV_PATH)
	MESSAGE(FATAL_ERROR "gcov not found! Aborting...")
ENDIF() # NOT GCOV_PATH

IF(NOT CMAKE_COMPILER_IS_GNUCXX)
	# Clang version 3.0.0 and greater now supports gcov as well.
	MESSAGE(WARNING "Compiler is not GNU gcc! Clang Version 3.0.0 and greater supports gcov as well, but older versions don't.")
	
	IF(NOT "${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang")
		MESSAGE(FATAL_ERROR "Compiler is not GNU gcc! Aborting...")
	ENDIF()
ENDIF() # NOT CMAKE_COMPILER_IS_GNUCXX


# Param _targetname     The name of new the custom make target
# Param _outputname     lcov output is generated as _outputname.info
#                       HTML report is generated in _outputname/index.html
# Param _testrunner     The name of the target which runs the tests.
#						MUST return ZERO always, even on errors. 
#						If not, no coverage report will be created!
# Optional fourth parameter is passed as arguments to _testrunner
#   Pass them in list form, e.g.: "-j;2" for -j 2
FUNCTION(SETUP_TARGET_FOR_COVERAGE _targetname _outputname _testrunner)

	IF(NOT LCOV_PATH)
		MESSAGE(STATUS "lcov not found! Cannot generate HTML coverage")
	ENDIF() # NOT LCOV_PATH

	IF(NOT GENHTML_PATH)
		MESSAGE(STATUS "genhtml not found! Cannot generate HTML coverage")
	ENDIF() # NOT GENHTML_PATH

	# Setup target
	ADD_CUSTOM_TARGET(${_targetname}
		
		# Cleanup lcov
		${LCOV_PATH} --directory . --zerocounters
		
		# Run tests
		COMMAND ${_testrunner} ${ARGV3}
		
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

ENDFUNCTION() # SETUP_TARGET_FOR_COVERAGE

# Param _targetname     The name of new the custom make target
# Param _outputname     cobertura output is generated as _outputname (e.g. coverage.xml)
# Param _testrunner     The name of the target which runs the tests
# Optional fourth parameter is passed as arguments to _testrunner
#   Pass them in list form, e.g.: "-j;2" for -j 2
FUNCTION(SETUP_TARGET_FOR_COVERAGE_COBERTURA _targetname _outputname _testrunner)

	# IF(NOT PYTHON_EXECUTABLE)
	# 	MESSAGE(FATAL_ERROR "Python not found! Aborting...")
	# ENDIF() # NOT PYTHON_EXECUTABLE

	IF(NOT GCOVR_PATH)
		MESSAGE(FATAL_ERROR "gcovr not found! Aborting...")
	ENDIF() # NOT GCOVR_PATH

	ADD_CUSTOM_TARGET(${_targetname}

		# Run tests
		${_testrunner} ${ARGV3}

		# Running gcovr
		COMMAND ${GCOVR_PATH} -x --xml-pretty -r ${CMAKE_SOURCE_DIR} -e '${CMAKE_SOURCE_DIR}/tests/'  -o ${_outputname}
		WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
		COMMENT "Running gcovr to produce Cobertura code coverage report."
	)

	# Show info where to find the report
	ADD_CUSTOM_COMMAND(TARGET ${_targetname} POST_BUILD
		COMMAND ;
		COMMENT "Cobertura code coverage report saved in ${_outputname}."
	)

ENDFUNCTION() # SETUP_TARGET_FOR_COVERAGE_COBERTURA




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
