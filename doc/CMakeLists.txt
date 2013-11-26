
find_package(Doxygen)

if(DOXYGEN_FOUND)
  message("GIT_DESCRIBE: ${GIT_DESCRIBE}")
  configure_file(${CMAKE_CURRENT_SOURCE_DIR}/doxy/flow123d.doxy.in 
    ${CMAKE_CURRENT_BINARY_DIR}/online-doc/Doxyfile 
  )
  add_custom_target(doc-doxy
    COMMAND ${DOXYGEN_EXECUTABLE} ${CMAKE_CURRENT_BINARY_DIR}/online-doc/Doxyfile
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/online-doc
    COMMENT "Generating API documentation with Doxygen" VERBATIM
  )
  set_directory_properties(PROPERTIES ADDITIONAL_MAKE_CLEAN_FILES "${CMAKE_BINARY_DIR}/doc/online-doc/flow123d;${CMAKE_BINARY_DIR}/doc/online-doc/flow123d.tag")
endif(DOXYGEN_FOUND)
