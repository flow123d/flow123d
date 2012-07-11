
set(CMAKE_INCLUDE_CURRENT_DIR ON)

set(TPeg_SOURCES 
    f123dTPeg.cc 
    
    yaml/yamlReader.cc 
    
    mesh/Physical.cc 
    mesh/Node.cc 
    mesh/Element.cc 
    mesh/EPoint.cc 
    mesh/ELine.cc 
    mesh/ETriangle.cc 
    mesh/ETetrahedron.cc 
    
    boundary/BoundaryData.cc 
    boundary/Boundary.cc)


if (YAML_FOUND)
   include_directories(${YAML_INCLUDE_DIR})
   add_executable(tpeg ${TPeg_SOURCES} )
   set_target_properties(tpeg PROPERTIES 
      RUNTIME_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/bin
   )
   target_link_libraries(tpeg ${YAML_LIBRARY})
endif (YAML_FOUND)

