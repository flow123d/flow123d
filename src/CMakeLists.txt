
set(CMAKE_INCLUDE_CURRENT_DIR ON)
message(STATUS "Include directories: ${Flow123d_SOURCE_DIR}")

set(libFlow123d_SRCS
    flow/darcy_flow_mh_output.cc
    flow/darcy_flow_mh.cc
    flow/local_matrix.cc
    mesh/boundaries.cc
    mesh/edges.cc
    mesh/elements.cc
    mesh/ini_constants_mesh.cc
    mesh/mesh.cc
    mesh/neighbours.cc
    mesh/nodes.cc
    mesh/sides.cc
    mesh/topology.cc
    reaction/linear_reaction.cc
    semchem/semchem_interface.cc
    system/math_fce.cc
    system/par_distribution.cc
    system/sys_function_stack.cpp
    system/sys_profiler.cc
    system/system.cc
    constantdb.cc
    convert.cc
    la_linsys.cc
    la_schur.cc
#    main.cc
    materials.cc
    msh_gmshreader.cc
    output.cc
    problem.cc
    read_ini.cc    
    solve.cc
    sparse_graph.cc
    time_governor.cc
    transport.cc
    xio.cc
)

set(libFlow123d_C_SRCS
    semchem/che_semchem.c
    semchem/cti_ichnew.c
)

foreach(file ${libFlow123d_C_SRCS})
  set(libFlow123d_SRCS ${libFlow123d_SRCS} ${file}c)
  configure_file(${file} ${file}c COPYONLY)
endforeach(file)

# add build targets
link_directories(${Flow123dLinkDirectories})
add_library(Flow123d ${libFlow123d_SRCS})
set_target_properties(Flow123d PROPERTIES
  ARCHIVE_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/lib
)

add_executable(flow123d main.cc)
add_custom_command(TARGET flow123d POST_BUILD 
  COMMAND ln -s ${CMAKE_SOURCE_DIR}/bin/${CURRENT_FLOW} ${CMAKE_BINARY_DIR}/bin/current_flow
)
set_directory_properties(PROPERTIES ADDITIONAL_MAKE_CLEAN_FILES ${CMAKE_BINARY_DIR}/bin/current_flow)
target_link_libraries(flow123d ${PETSC_LIBRARIES} ${PARMETIS_LIBRARIES} Flow123d)
set_target_properties(flow123d PROPERTIES
   RUNTIME_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/bin
   OBJECT_DIR ${PROJECT_BINARY_DIR}/build/obj
)