
enable_testing()

# Set timeout for one unit test
set(TEST_TIMEOUT 30) # 30 s
# Set timeout for one unit test running under valgrind
set(VALGRIND_TIMEOUT 300) # 5 min

# build gtest itself
set(GTEST_SOURCE_DIR "${CMAKE_SOURCE_DIR}/third_party/gtest-1.7.0")
set(GTEST_BINARY_DIR "${CMAKE_BINARY_DIR}/third_party/gtest-1.7.0")
include_directories( ${GTEST_SOURCE_DIR}/include ${CMAKE_SOURCE_DIR}/unit_tests )
add_subdirectory( ${GTEST_SOURCE_DIR} ${GTEST_BINARY_DIR} )
message(STATUS "Gtest tests allowed")

# sources for our mpi gtest extension
set(mpi_gtest_sources  ${CMAKE_SOURCE_DIR}/unit_tests/_gtest_extensions/gtest_mpi.cc)
# object library used by tests, but can be build once
add_library(gtest_mpi_obj OBJECT ${mpi_gtest_sources})
add_dependencies(gtest_mpi_obj gtest gtest_main)

# Defines macro for absolute path to unit_tests, so that unit tests can read
# some test files in its source directory, while the binary of the test is created
# elsewhere.
add_definitions(-DUNIT_TESTS_SRC_DIR="${CMAKE_SOURCE_DIR}/unit_tests")

# target to make all unit tests
add_custom_target(all-unit-test)

######################################################################################################
# Macro to make default targets and operations in current unit tests subdirectory, namely:
# - create sym link to the generated Makefile
# - make all-test and all-bench targets
# - add required libraries
# Should be called from subdirectories
macro(add_test_directory libs)
  
  string(REPLACE "${CMAKE_SOURCE_DIR}/unit_tests/" "" subdir "${CMAKE_CURRENT_SOURCE_DIR}")
  #message("IN: ${CMAKE_CURRENT_SOURCE_DIR} SUBDIR: ${subdir}")
  set(super_test "${subdir}-test")
  add_custom_target( ${super_test} )
  set(super_test_libs ${libs})
  #message("SUBDIR LIBS: ${libs}")
  
  if (ALLOW_DUPLICATE_CUSTOM_TARGETS)
    #message("allow common targets")
    
    # default all_test target
    add_custom_target(all-test
                    DEPENDS ${super_test}
                    )
  
    # default all_bench target
    add_custom_target(all-bench
                    #DEPENDS ${test_name}
                    )
  endif()
  
  add_dependencies(all-unit-test ${super_test})
  
  make_link_in_source("${CMAKE_CURRENT_SOURCE_DIR}/Makefile")    
endmacro()



#################################
# MACRO define_test
#
# Creates build targets for particular source file containing one test.
#
# For given parameter 'class_name', it assumes class_name_test.cpp source
# and creates:
# target 'class_name_test' which builds corresponding binary and run the test
# target 'class_name_valgrind' which builds the binary and runs it under valgrind
# use libraries in ${super_test_libs}
#
# We want to run tests from source dir, but since there is no portable way 
# how to make link from source dir to Makefile generated in build dir 
# I let main makefile to make such links.
macro(define_test class_name)
  set(test_source "${class_name}_test.cpp")
  set(test_binary "${class_name}_test_bin")
  set(test_name "${class_name}-test")
  set(test_valgrind "${class_name}-valgrind")
  #message(STATUS "adding test: ${test_binary} cxx flags: ${CMAKE_CXX_FLAGS}")
  
  # executable rule
  add_executable(${test_binary} EXCLUDE_FROM_ALL  ${test_source} )
  #set_target_properties(${test_binary} PROPERTIES COMPILE_FLAGS "${CMAKE_CXX_FLAGS}")
  target_link_libraries(${test_binary} ${super_test_libs} gtest gtest_main)
  add_dependencies(${test_binary} update_source_tree)


  
  # rule to run the test
  add_custom_target(${test_name} 
      COMMAND ${CMAKE_SOURCE_DIR}/bin/time_limit.sh -t ${TEST_TIMEOUT} ${CMAKE_CURRENT_BINARY_DIR}/${test_binary} --gtest_output="xml:${CMAKE_CURRENT_BINARY_DIR}/${test_name}.xml"
      DEPENDS ${test_binary}
      WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})

  # rule to run the test under valgrind
  add_custom_target(${test_valgrind} 
      COMMAND ${CMAKE_SOURCE_DIR}/bin/time_limit.sh -t ${VALGRIND_TIMEOUT} valgrind ${CMAKE_CURRENT_BINARY_DIR}/${test_binary} --gtest_output="xml:${CMAKE_CURRENT_BINARY_DIR}/${test_valgrind}.xml"
      DEPENDS ${test_binary}
      WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})
      
  # register the test into super_test
  add_dependencies(${super_test} ${test_name})

  # let ctest call make with appropriate  target
  add_test(NAME ${test_name} COMMAND ${CMAKE_MAKE_PROGRAM} ${test_name})
endmacro()

#################################
# MACRO define_mpi_test
#
# Creates build targets for particular source file containing one test which should run in parallel.
#
# It has one more parameter for number of processors N. It provides similar targets as the macro 'define_test'
# like: 'class_name_N_test' and 'class_name_N_valgrind' which builds the test binary and run it in parallel using bin/mpiexec
# 
# To get meaningful output from the parallel tests you should include <gtest_flow.hh> instead of <gtest/gtest.h>
macro(define_mpi_test class_name n_proc)
  set(test_source "${class_name}_test.cpp")
  set(test_binary "${class_name}_test_bin")
  set(test_name "${class_name}-${n_proc}-test")
  set(test_valgrind "${class_name}-${n_proc}-valgrind")
  #message(STATUS "adding mpi test: ${test_binary} on nprocs: ${n_proc} cxx flags: ${CMAKE_CXX_FLAGS}")
  
  # executable rule
  if(TARGET ${test_binary})
  else() 
  	add_executable(${test_binary} EXCLUDE_FROM_ALL   ${test_source} $<TARGET_OBJECTS:gtest_mpi_obj>) 
  	target_link_libraries(${test_binary} ${super_test_libs} gtest gtest_main  ${MPI_CXX_LIBRARIES})
  	add_dependencies(${test_binary} mpiexec_link)      # make
  endif()	

  # rule to run the test
  # set(RUN_COMMAND ulimit -t 300\; )
  # message(STATUS "BASH_COMMAND=${BASH_COMMAND}")  
  add_custom_target(${test_name}
      COMMAND ${CMAKE_SOURCE_DIR}/bin/time_limit.sh -t ${TEST_TIMEOUT} ${CMAKE_BINARY_DIR}/bin/mpiexec -np ${n_proc} ${CMAKE_CURRENT_BINARY_DIR}/${test_binary} --gtest_output="xml:${CMAKE_CURRENT_BINARY_DIR}/${test_name}.xml"
      DEPENDS ${test_binary}
      WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})
      	      

  # rule to run the test under valgrind
  add_custom_target(${test_valgrind} 
      COMMAND ${CMAKE_SOURCE_DIR}/bin/time_limit.sh -t ${VALGRIND_TIMEOUT} ${CMAKE_BINARY_DIR}/bin/mpiexec -np ${n_proc} valgrind ${CMAKE_CURRENT_BINARY_DIR}/${test_binary} --gtest_output="xml:${CMAKE_CURRENT_BINARY_DIR}/${test_valgrind}.xml"
      DEPENDS ${test_binary}
      WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})
      
  # register the test into super_test
  add_dependencies(${super_test} ${test_name})

  # let ctest call make with appropriate  target
  add_test(NAME ${test_name} COMMAND ${CMAKE_MAKE_PROGRAM} ${test_name})
endmacro()



#######################################################################################################
# add subdirectories
# - has to be done manually to get correct dependencies on CMake's list files.


add_subdirectory("scripts")
add_subdirectory("input")
add_subdirectory("interpolation")
add_subdirectory("la")
add_subdirectory("system")
add_subdirectory("fields")
add_subdirectory("fem")
add_subdirectory("mesh")
add_subdirectory("coupling")
add_subdirectory("output")


#################################################################################################
# Targets for Jenkis Jobs Builder

# Make list of integration test names from list of subdirectories in the test directory. 
# This should be manually finished before it can be used by use it in 
add_custom_target(test-list
                  COMMAND ${CMAKE_COMMAND} -DDIRECTORY="${CMAKE_SOURCE_DIR}/tests" -DOUTPUT="${CMAKE_CURRENT_SOURCE_DIR}/test_list" -P ${CMAKE_CURRENT_SOURCE_DIR}/list_subdirs.cmake
                  )


# Following two targets assume installation of Jenkins jobs Builder.
# https://github.com/openstack-infra/jenkins-job-buildera
#
# For installation follow README.md.

# Rule to create XML definitions of Jenkins jobs and upload them to the server.
set(JENKINS_JOBS_FILE ${CMAKE_CURRENT_BINARY_DIR}/jenkins-jobs-config.yaml)
add_custom_target(jenkins-upload
                  COMMAND jenkins-jobs --ignore-cache update ${JENKINS_JOBS_FILE} 
                  DEPENDS ${JENKINS_JOBS_FILE}
                  WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
                  )

# Rule to create XML definitions of Jenkins jobs in local directory. (for testing)
add_custom_target(jenkins-test
                  COMMAND jenkins-jobs test ${JENKINS_JOBS_FILE} -o ${CMAKE_CURRENT_BINARY_DIR}                 
                  DEPENDS ${JENKINS_JOBS_FILE}
                  WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
                  )

# Rule how to make ${JENKINS_JOBS_FILE} from ./jenkins_jobs_template.yaml
add_custom_command(OUTPUT ${JENKINS_JOBS_FILE}
                  COMMAND ${CMAKE_COMMAND} -DINPUT=jenkins_jobs_template.yaml -DOUTPUT=${JENKINS_JOBS_FILE} -P ${CMAKE_CURRENT_SOURCE_DIR}/expand_placeholders.cmake
                  DEPENDS jenkins_jobs_template.yaml
                  WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
                  )

