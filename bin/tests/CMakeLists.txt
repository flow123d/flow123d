  # target to compile tester_mpi.cc
  add_executable(tester_mpi EXCLUDE_FROM_ALL tester_mpi.cc) 
  target_link_libraries(tester_mpi ${MPI_CXX_LIBRARIES})
  add_dependencies(tester_mpi mpiexec_link)      # need mpiexec link compatible with MPI lib

  message(STATUS "flow123d.sh test: ${CMAKE_CURRENT_BINARY_DIR}")
  # copy flow123d.sh and its test to testing build dir
  #file(COPY ../config
  #     DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/test_dir)     

  # target to call flow123d.sh.test
  add_custom_target(test_flow123d_sh
      COMMAND ${CMAKE_COMMAND} -E make_directory ${CMAKE_CURRENT_BINARY_DIR}/test_dir
      COMMAND ${CMAKE_COMMAND} -E copy ${CMAKE_SOURCE_DIR}/bin/time_limit.sh ${CMAKE_CURRENT_BINARY_DIR}/test_dir
      COMMAND ${CMAKE_COMMAND} -E copy ${CMAKE_SOURCE_DIR}/bin/flow123d.sh ${CMAKE_CURRENT_BINARY_DIR}/test_dir
      COMMAND ${CMAKE_COMMAND} -E copy ${CMAKE_BINARY_DIR}/mpiexec ${CMAKE_CURRENT_BINARY_DIR}/test_dir
      COMMAND ${CMAKE_COMMAND} -E copy_directory ${CMAKE_SOURCE_DIR}/bin/config ${CMAKE_CURRENT_BINARY_DIR}/test_dir
      COMMAND ${CMAKE_SOURCE_DIR}/bin/tests/flow123d.sh.test 
      DEPENDS tester_mpi
      WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})
