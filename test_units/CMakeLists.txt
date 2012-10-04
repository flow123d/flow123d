
enable_testing()

# build gtest itself
SET(GTEST_ROOT "${CMAKE_SOURCE_DIR}/third_party/gtest-1.6.0")
include_directories( ${GTEST_ROOT}/include ${CMAKE_SOURCE_DIR}/test_units)
add_subdirectory( ${GTEST_ROOT} ${GTEST_ROOT})
message(STATUS "Gtest tests allowed")

set(mpi_gtest_sources  ${CMAKE_SOURCE_DIR}/test_units/gtest_mpi.cc)
set(CMAKE_CXX_FLAGS "-g") 

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
  set(test_name "${class_name}_test")
  set(test_valgrind "${class_name}_valgrind")
  
  # executable rule
  add_executable(${test_binary} EXCLUDE_FROM_ALL  ${test_source} ) 
  target_link_libraries(${test_binary} ${super_test_libs} gtest gtest_main)

  # rule to run the test
  add_custom_target(${test_name} 
      COMMAND ${test_binary}
      DEPENDS ${test_binary}
      WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})

  # rule to run the test under valgrind
  add_custom_target(${test_valgrind} 
      COMMAND valgrind ${CMAKE_CURRENT_BINARY_DIR}/${test_binary}
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
  set(test_name "${class_name}_${n_proc}_test")
  set(test_valgrind "${class_name}_${n_proc}_valgrind")
  
  # executable rule
  if(TARGET ${test_binary})
  else() 
  	add_executable(${test_binary} EXCLUDE_FROM_ALL   ${test_source} ${mpi_gtest_sources}) 
  	target_link_libraries(${test_binary} ${super_test_libs} gtest gtest_main  ${MPI_CXX_LIBRARIES})
  	add_dependencies(${test_binary} mpiexec_link)      # make
  endif()	

  # rule to run the test
  add_custom_target(${test_name} 
      COMMAND ${CMAKE_BINARY_DIR}/mpiexec -np ${n_proc} ${CMAKE_CURRENT_BINARY_DIR}/${test_binary} --gtest_color=no 
      DEPENDS ${test_binary}
      WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})
      	      
      	      #| ${CMAKE_SOURCE_DIR}/bin/tests/filter_output.sh ${CMAKE_CURRENT_SOURCE_DIR}/mpi_out 
      	      

  # rule to run the test under valgrind
  add_custom_target(${test_valgrind} 
      COMMAND ${CMAKE_BINARY_DIR}/mpiexec -np ${n_proc} valgrind ${CMAKE_CURRENT_BINARY_DIR}/${test_binary} --gtest_color=yes
      DEPENDS ${test_binary}
      WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})
      
  # register the test into super_test
  add_dependencies(${super_test} ${test_name})

  # let ctest call make with appropriate  target
  add_test(NAME ${test_name} COMMAND ${CMAKE_MAKE_PROGRAM} ${test_name})
endmacro()

##############################
# test subdirs
add_subdirectory("input")
add_subdirectory("interpolation")
add_subdirectory("la")
add_subdirectory("system")



#add_executable(ioNameHandlerTest system/io_name_handler.cc)
#set_target_properties(ioNameHandlerTest PROPERTIES
#     RUNTIME_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/test_units/bin
#  )
#  target_link_libraries(ioNameHandlerTest ${PETSC_LIBRARIES} ${PARMETIS_LIBRARIES} ${GTEST_LIBRARIES} Flow123d)

#add_test(NAME io_name_handler 
#    COMMAND ioNameHandlerTest -o /output -i input)

