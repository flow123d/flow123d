
enable_testing()

# Set timeout for one unit test

enable_testing()

# Set timeout for one unit test
set(TEST_TIMEOUT 15) # [s] Hard limit for Jenkins. Unit tests should run under 10 s.
# Set timeout for one unit test running under valgrind
set(VALGRIND_TIMEOUT 15000) # 1000x

# build gtest itself
set(GTEST_SOURCE_DIR "${CMAKE_SOURCE_DIR}/third_party/gtest-1.10.0")
set(GTEST_BINARY_DIR "${CMAKE_BINARY_DIR}/third_party/gtest-1.10.0")
include_directories( ${GTEST_SOURCE_DIR}/googletest/include ${CMAKE_SOURCE_DIR}/unit_tests )
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
add_definitions(-DUNIT_TESTS_BIN_DIR="${CMAKE_BINARY_DIR}/unit_tests")

# target to make all unit tests
#add_custom_target(all-unit-test)

######################################################################################################
# Macro to make default targets and operations in current unit tests subdirectory, namely:
# - create sym link to the generated Makefile
# - make all-test and all-bench targets
# - add required libraries
# Should be called from subdirectories
macro(add_test_directory libs)

  string(REPLACE "${CMAKE_SOURCE_DIR}/unit_tests/" "" subdir "${CMAKE_CURRENT_SOURCE_DIR}")
  #message("IN: ${CMAKE_CURRENT_SOURCE_DIR} SUBDIR: ${subdir}")
  set(dir_test_target "${subdir}-test")
  set(dir_bench_target "${subdir}-bench")
  add_custom_target( ${dir_test_target} )
  add_custom_target( ${dir_bench_target} )
  set(dir_test_target_libs ${libs})
  #message("SUBDIR LIBS: ${libs}")

  if (ALLOW_DUPLICATE_CUSTOM_TARGETS)
    #message("allow common targets")

    # default all_test target
    add_custom_target(all-test
                    DEPENDS ${dir_test_target}
                    )

    # default all_bench target
    add_custom_target(all-bench
                    DEPENDS ${dir_bench_target}
                    )
  endif()

  #add_dependencies(all-unit-test ${dir_test_target})

  make_link_in_source("${CMAKE_CURRENT_SOURCE_DIR}/Makefile")
endmacro()



#################################
# MACRO define_test
# For given parameter 'class_name', it assumes class_name_test.cpp source
# and creates:
# target 'class_name_test' which builds corresponding binary and run the test
# target 'class_name_valgrind' which builds the binary and runs it under valgrind
#
# Creates build targets for particular source file containing one test.
#
# use libraries in ${dir_test_target_libs}
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

  set (extra_args ${ARGN})
  # Did we get any optional args?
  list(LENGTH extra_args extra_count)
  if (${extra_count} GREATER 0)
      list(GET extra_args 0 test_timeout)
  else () 
      set(test_timeout ${TEST_TIMEOUT})
  endif ()
  math(EXPR valgrind_timeout "${test_timeout} * 1000")
  
  # executable rule
  add_executable(${test_binary} EXCLUDE_FROM_ALL  ${test_source} )
  set_target_properties(${test_binary} PROPERTIES COMPILE_FLAGS "${CMAKE_CXX_FLAGS} -Wno-error")
  target_link_libraries(${test_binary} ${dir_test_target_libs} gtest gtest_main)
  add_dependencies(${test_binary} update_source_tree)



  # rule to run the test
  add_custom_target(${test_name}
      COMMAND ${CMAKE_SOURCE_DIR}/bin/time_limit.sh -t ${test_timeout} ${CMAKE_CURRENT_BINARY_DIR}/${test_binary} --gtest_output="xml:${CMAKE_CURRENT_BINARY_DIR}/${test_name}.xml"
      DEPENDS ${test_binary}
      WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})

  # rule to run the test under valgrind
  add_custom_target(${test_valgrind}
      COMMAND ${CMAKE_SOURCE_DIR}/bin/time_limit.sh -t ${valgrind_timeout} valgrind --num-callers=50 --suppressions=${CMAKE_SOURCE_DIR}/bin/python.supp ${CMAKE_CURRENT_BINARY_DIR}/${test_binary} --gtest_output="xml:${CMAKE_CURRENT_BINARY_DIR}/${test_valgrind}.xml"
      DEPENDS ${test_binary}
      WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})

  # register the test into dir_test_target
  add_dependencies(${dir_test_target} ${test_name})

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

  set (extra_args ${ARGN})
  # Did we get any optional args?
  list(LENGTH extra_args extra_count)
  if (${extra_count} GREATER 0)
      list(GET extra_args 0 test_timeout)
  else () 
      set(test_timeout ${TEST_TIMEOUT})
  endif ()
  math(EXPR valgrind_timeout "${test_timeout} * 1000")

  # executable rule
  if(TARGET ${test_binary})
  else()
  	add_executable(${test_binary} EXCLUDE_FROM_ALL   ${test_source} $<TARGET_OBJECTS:gtest_mpi_obj>)
  	target_link_libraries(${test_binary} ${dir_test_target_libs} gtest gtest_main  ${MPI_CXX_LIBRARIES})
  	# add_dependencies(${test_binary} mpiexec_link)      # make
  endif()
  set_target_properties(${test_binary} PROPERTIES COMPILE_FLAGS "${CMAKE_CXX_FLAGS} -Wno-error")

  # rule to run the test
  # set(RUN_COMMAND ulimit -t 300\; )
  # message(STATUS "BASH_COMMAND=${BASH_COMMAND}")
  add_custom_target(${test_name}
      COMMAND ${CMAKE_SOURCE_DIR}/bin/time_limit.sh -t ${test_timeout} ${CMAKE_BINARY_DIR}/bin/mpiexec -np ${n_proc} ${CMAKE_CURRENT_BINARY_DIR}/${test_binary} --gtest_output="xml:${CMAKE_CURRENT_BINARY_DIR}/${test_name}.xml"
      DEPENDS ${test_binary}
      WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})


  # rule to run the test under valgrind
  add_custom_target(${test_valgrind}
      COMMAND ${CMAKE_SOURCE_DIR}/bin/time_limit.sh -t ${valgrind_timeout} ${CMAKE_BINARY_DIR}/bin/mpiexec -np ${n_proc} valgrind --num-callers=50 --suppressions=${CMAKE_SOURCE_DIR}/bin/python.supp ${CMAKE_CURRENT_BINARY_DIR}/${test_binary} --gtest_output="xml:${CMAKE_CURRENT_BINARY_DIR}/${test_valgrind}.xml"
      DEPENDS ${test_binary}
      WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})

  # register the test into dir_test_target
  add_dependencies(${dir_test_target} ${test_name})

  # let ctest call make with appropriate  target
  add_test(NAME ${test_name} COMMAND ${CMAKE_MAKE_PROGRAM} ${test_name})
endmacro()

#################################
# MACRO define_mpi_benchmark
#
# Like define_mpi_test, but with following differences:
# - is executed in the directory of the benchmark binary
# - is not part of all_test
# - accepts additional parameter specifying the python postprocessing script for collecting and possibly publissing the benchmark results
#   namely from the profiler output
#
# It has one more parameter for number of processors N. It provides similar targets as the macro 'define_test'
# like: 'class_name_N_test' and 'class_name_N_valgrind' which builds the test binary and run it in parallel using bin/mpiexec
#
# To get meaningful output from the parallel tests you should include <gtest_flow.hh> instead of <gtest/gtest.h>
macro(define_mpi_benchmark class_name n_proc postproc_script)
  set(test_source "${class_name}_bench.cpp")
  set(test_binary "${class_name}_bench_bin")
  set(test_name "${class_name}-${n_proc}-bench")
  set(test_valgrind "${class_name}-${n_proc}-valgrind")
  #message(STATUS "adding mpi test: ${test_binary} on nprocs: ${n_proc} cxx flags: ${CMAKE_CXX_FLAGS}")

  set (extra_args ${ARGN})
  # Did we get any optional args?
  list(LENGTH extra_args extra_count)
  if (${extra_count} GREATER 0)
      list(GET extra_args 0 test_timeout)
  else () 
      set(test_timeout ${TEST_TIMEOUT})
  endif ()
  math(EXPR valgrind_timeout "${test_timeout} * 1000")

  # executable rule
  if(TARGET ${test_binary})
  else()
  	add_executable(${test_binary} EXCLUDE_FROM_ALL   ${test_source} $<TARGET_OBJECTS:gtest_mpi_obj>)
  	target_link_libraries(${test_binary} ${dir_test_target_libs} gtest gtest_main  ${MPI_CXX_LIBRARIES})
  	# add_dependencies(${test_binary} mpiexec_link)      # make
  endif()
  set_target_properties(${test_binary} PROPERTIES COMPILE_FLAGS "${CMAKE_CXX_FLAGS} -Wno-error")

  # rule to run the test
  # set(RUN_COMMAND ulimit -t 300\; )
  # message(STATUS "BASH_COMMAND=${BASH_COMMAND}")
  add_custom_target(${test_name}
      # BENCHMARK_N_RUNS variable can be used to specify number of repetition
      COMMAND ${CMAKE_SOURCE_DIR}/bin/run_benchmark.sh -t=${test_timeout} -np=${n_proc}  ${CMAKE_CURRENT_BINARY_DIR}/${test_binary} 
      DEPENDS ${test_binary}
      WORKING_DIRECTORY ${CMAKE_SOURCE_DIR})


  # rule to run the test under valgrind
  #add_custom_target(${test_valgrind}
  #    COMMAND ${CMAKE_SOURCE_DIR}/bin/time_limit.sh -t ${valgrind_timeout} ${CMAKE_BINARY_DIR}/bin/mpiexec -np ${n_proc} valgrind --num-callers=50 --suppressions=${CMAKE_SOURCE_DIR}/bin/#python.supp ${CMAKE_CURRENT_BINARY_DIR}/${test_binary} --gtest_output="xml:${CMAKE_CURRENT_BINARY_DIR}/${test_valgrind}.xml"
  #    DEPENDS ${test_binary}
  #    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})

  # register the test into dir_test_target
  add_dependencies(${dir_bench_target} ${test_name})

  # let ctest call make with appropriate  target
  add_test(NAME ${test_name} COMMAND ${CMAKE_MAKE_PROGRAM} ${test_name})
endmacro()


#######################################################################################################
# add subdirectories
# - has to be done manually to get correct dependencies on CMake's list files.


add_subdirectory("test_scripts")
add_subdirectory("scripts")
add_subdirectory("input")
add_subdirectory("flow")
add_subdirectory("tools")
add_subdirectory("la")
add_subdirectory("system")
add_subdirectory("fields")
add_subdirectory("fem")
add_subdirectory("mesh")
add_subdirectory("intersection")
add_subdirectory("coupling")
add_subdirectory("output")
add_subdirectory("dealii")


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
                  COMMAND jenkins-jobs --ignore-cache --conf ${CMAKE_SOURCE_DIR}/config/jenkins.ini update ${JENKINS_JOBS_FILE}
                  DEPENDS ${JENKINS_JOBS_FILE}
                  WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
                  )

# Rule to create XML definitions of Jenkins jobs in local directory. (for testing)
add_custom_target(jenkins-test
                  COMMAND jenkins-jobs --conf ${CMAKE_SOURCE_DIR}/config/jenkins.ini test ${JENKINS_JOBS_FILE} -o ${CMAKE_CURRENT_BINARY_DIR}
                  DEPENDS ${JENKINS_JOBS_FILE}
                  WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
                  )

# Rule how to make ${JENKINS_JOBS_FILE} from ./jenkins_jobs_template.yaml
add_custom_command(OUTPUT ${JENKINS_JOBS_FILE}
                  COMMAND ${CMAKE_COMMAND} -DINPUT=jenkins_jobs_template.yaml -DOUTPUT=${JENKINS_JOBS_FILE} -P ${CMAKE_CURRENT_SOURCE_DIR}/expand_placeholders.cmake
                  DEPENDS jenkins_jobs_template.yaml
                  WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
                  )
