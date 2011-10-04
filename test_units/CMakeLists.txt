
if(GTEST_FOUND)
  message(STATUS "Gtest tests allowed")
  add_executable(ioNameHandlerTest system/io_name_handler.cc)

  set_target_properties(ioNameHandlerTest PROPERTIES
     RUNTIME_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/test_units/bin
  )
  target_link_libraries(ioNameHandlerTest ${PETSC_LIBRARIES} ${PARMETIS_LIBRARIES} ${GTEST_LIBRARIES} Flow123d)

  add_test(NAME io_name_handler 
    COMMAND ioNameHandlerTest -o /output -i input)
endif(GTEST_FOUND)

if(CppUnit_FOUND)
  message(STATUS "CppUnit tests allowed")


endif(CppUnit_FOUND)