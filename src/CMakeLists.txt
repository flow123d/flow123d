
set(CMAKE_INCLUDE_CURRENT_DIR ON)
message(STATUS "Include directories: ${Flow123d_SOURCE_DIR}")

set(libFlow123d_SRCS
    boundaries.cc
#    btc.cc
#    calculation.cc
#    concentrations.cc
    constantdb.cc
    convert.cc
    darcy_flow_mh_output.cc
    darcy_flow_mh.cc
    edges.cc
    elements.cc
  #  hashes.cc
    la_linsys.cc
    la_schur.cc
    local_matrix.cc
    main.cc
    materials.cc
    math_fce.cc
    mesh.cc
    msh_gmshreader.cc
    neighbours.cc
    nodes.cc
    output.cc
    par_distribution.cc
#    ppfcs.cc
    problem.cc
#    reaction.cc
    read_ini.cc
    sides.cc
    solve.cc
    sparse_graph.cc
    sys_function_stack.cpp
    sys_profiler.cc
    system.cc
    time_governor.cc
    topology.cc
#    transport_bcd.cc
    transport.cc
    xio.cc
    reaction/linear_reaction.cc
    mesh/ini_constants_mesh.cc
    semchem/semchem_interface.cc
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
add_executable(flow123d ${libFlow123d_SRCS})
target_link_libraries(flow123d ${PETSC_LIBRARIES} ${PARMETIS_LIBRARIES})
set_target_properties(flow123d PROPERTIES
   RUNTIME_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/bin
   OBJECT_DIR ${PROJECT_BINARY_DIR}/build/obj
)