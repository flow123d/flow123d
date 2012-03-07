
enable_testing()

# build gtest itself
SET(GTEST_ROOT "${CMAKE_SOURCE_DIR}/third_party/gtest-1.6.0")
include_directories( ${GTEST_ROOT}/include )
add_subdirectory( ${GTEST_ROOT} ${GTEST_ROOT})


message(STATUS "Gtest tests allowed")

add_subdirectory(input)

#add_executable(ioNameHandlerTest system/io_name_handler.cc)
#set_target_properties(ioNameHandlerTest PROPERTIES
#     RUNTIME_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/test_units/bin
#  )
#  target_link_libraries(ioNameHandlerTest ${PETSC_LIBRARIES} ${PARMETIS_LIBRARIES} ${GTEST_LIBRARIES} Flow123d)

#add_test(NAME io_name_handler 
#    COMMAND ioNameHandlerTest -o /output -i input)

