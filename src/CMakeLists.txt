# 
# Copyright (C) 2007 Technical University of Liberec.  All rights reserved.
#
# Please make a following refer to Flow123d on your project site if you use the program for any purpose,
# especially for academic research:
# Flow123d, Research Centre: Advanced Remedial Technologies, Technical University of Liberec, Czech Republic
#
# This program is free software; you can redistribute it and/or modify it under the terms
# of the GNU General Public License version 3 as published by the Free Software Foundation.
# 
# This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; 
# without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
# See the GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License along with this program; if not,
# write to the Free Software Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA 021110-1307, USA.
#
# $Id$
# $Revision$
# $LastChangedBy$
# $LastChangedDate$
#


# when the precompiled headers are used in the future 
#if (${PRECOMPILED_HEADERS})
#    include(CMake/Modules/PrecompiledHeaders.cmake)
#endif()

add_subdirectory("input/json_spirit")

set(CMAKE_INCLUDE_CURRENT_DIR ON)
#message(STATUS "Include directories: ${flow123d_includes}")
#include_directories(${flow123d_includes})




# select source files that compose whole flow123d library
# todo: split into parts, e.g. system, output, input, mesh, darcy_mh, transport, .. 
# every library then can be testetd independently
set(sources_lib_system
    system/python_loader.cc
    system/math_fce.cc
    system/sys_function_stack.cpp
    system/sys_profiler.cc
    system/system.cc
    system/xio.cc
    system/exceptions.cc
    io/read_ini.cc
    system/file_path.cc
    system/tokenizer.cc
)

set(sources_lib_input
  
    input/type_base.cc
    input/type_record.cc
    input/type_selection.cc
    input/storage.cc
    input/json_to_storage.cc
    input/type_output.cc
    
    input/accessors.cc
)


set(sources_lib_la
	la/distribution.cc
	la/local_to_global_map.cc
    la/schur.cc
    la/solve.cc
    la/linsys.cc
    la/linsys_BDDC.cc
    la/linsys_PETSC.cc
    la/sparse_graph.cc
)

set(sources_mesh_lib
    mesh/region.cc
    mesh/msh_gmshreader.cc
    mesh/boundaries.cc
    mesh/edges.cc
    mesh/elements.cc
    mesh/mesh.cc
    mesh/neighbours.cc
    mesh/nodes.cc
    mesh/sides.cc
    mesh/intersection.cc
    mesh/ref_element.cc
    mesh/partitioning.cc
    
    mesh/bounding_box.cc
    mesh/bih_node.cc
    mesh/bih_tree.cc
    
    mesh/ngh/src/abscissa.cpp
    mesh/ngh/src/bisector.cpp
    mesh/ngh/src/intersection.cpp
    mesh/ngh/src/intersectionLocal.cpp
    mesh/ngh/src/mathfce.cpp
    mesh/ngh/src/matrix.cpp
    mesh/ngh/src/plain.cpp
    mesh/ngh/src/point.cpp
    mesh/ngh/src/polygon.cpp
    mesh/ngh/src/system.cpp
    mesh/ngh/src/tetrahedron.cpp
    mesh/ngh/src/triangle.cpp
    mesh/ngh/src/vector.cpp
    mesh/ngh/src/vertex.cpp
)

set(sources_fem_lib
    quadrature/quadrature_lib.cc

    fem/fe_p.cc
    fem/dofhandler.cc
    fem/fe_values.cc
    fem/finite_element.cc
    fem/mapping_p1.cc
    
)





set(sources_lib_field
    fields/field_values.cc
    fields/field_base.cc
    fields/field_other.cc
)


set (sources_lib_coupling

    coupling/time_marks.cc
    coupling/time_governor.cc

   coupling/equation.cc
       flow/old_bcd.cc
   # coupling/hc_explicit_sequential.cc
)

set(sources_lib_non_system

    flow/mh_fe_values.cc
    flow/mh_dofhandler.cc
    flow/darcy_flow_mh_output.cc
    flow/darcy_flow_mh.cc


    io/output.cc
    io/output_vtk.cc
    io/output_msh.cc
    
    reaction/reaction.cc
    reaction/linear_reaction.cc
    reaction/pade_approximant.cc
    reaction/sorption.cc
    reaction/isotherm.cc
    semchem/semchem_interface.cc

    transport/transport.cc
    transport/transport_operator_splitting.cc
    transport/transport_dg.cc

    coupling/hc_explicit_sequential.cc
)



# make separate library for Semchem (God save us!)
set(sources_lib_semchem
    semchem/che_semchem.cc
    semchem/che_read.cc
#    semchem/cti_ichnew.c
)



# add individual libraries

if (mpiuni_sources) 
  add_library(mpiuni_lib ${mpiuni_sources})
endif()

# system_lib
add_library(system_lib ${sources_lib_system})
target_link_libraries(system_lib ${MPI_CXX_LIBRARIES} ${Boost_LIBRARIES} ${PYTHON_LIBRARIES})


# when the precompiled headers are used in the future
#if (${PRECOMPILED_HEADERS})
#    ADD_PRECOMPILED_HEADER(system_lib "pre_system.hh")
#endif()



# input_lib
add_library(input_lib ${sources_lib_input})
target_link_libraries(input_lib 
    system_lib json_spirit)

    
# la_lib
add_library(la_lib ${sources_lib_la})
target_link_libraries(la_lib 
    input_lib system_lib
    armadillo ${PETSC_LIBRARIES}  ${Boost_LIBRARIES} ${BDDCML_LIBRARIES})

# mesh_lib
add_library(mesh_lib ${sources_mesh_lib})
target_link_libraries(mesh_lib
    input_lib system_lib la_lib armadillo ${Boost_LIBRARIES})
    
# fem_lib    
add_library(fem_lib ${sources_fem_lib})
target_link_libraries(fem_lib
    system_lib mesh_lib armadillo ${Boost_LIBRARIES})

# field_lib  
add_library(field_lib ${sources_lib_field})
target_link_libraries(field_lib
    system_lib input_lib fparser mesh_lib armadillo  ${Boost_LIBRARIES})
 
# coupling_lib  
add_library(coupling_lib ${sources_lib_coupling})
target_link_libraries(coupling_lib
    system_lib input_lib field_lib mesh_lib)
  
# all remaining sources 
add_library(flow123d_lib ${sources_lib_non_system})
target_link_libraries(flow123d_lib 
    fem_lib mesh_lib la_lib semchem  input_lib system_lib
    coupling_lib fparser  armadillo ${PETSC_LIBRARIES}  ${Boost_LIBRARIES})
set_target_properties(flow123d_lib PROPERTIES
    ARCHIVE_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/lib
)

# add target for libsemchem
add_library(semchem ${sources_lib_semchem})
set_target_properties(semchem PROPERTIES
  ARCHIVE_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/lib
)







set_directory_properties(PROPERTIES ADDITIONAL_MAKE_CLEAN_FILES "${CMAKE_BINARY_DIR}/bin/current_flow;${CMAKE_BINARY_DIR}/bin/mpiexec")


