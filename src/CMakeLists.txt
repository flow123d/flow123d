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

set(CMAKE_INCLUDE_CURRENT_DIR ON)
#message(STATUS "Include directories: ${flow123d_includes}")
#include_directories(${flow123d_includes})

# select source files that compose whole flow123d library
# todo: split into parts, e.g. system, output, input, mesh, darcy_mh, transport, .. 
# every library then can be testetd independently
set(sources_lib_non_system
    fem/dofhandler.cc
    fem/fe_values.cc
    fem/finite_element.cc
    fem/mapping_p1.cc

    flow/local_matrix.cc
    flow/darcy_flow_mh_output.cc
    flow/darcy_flow_mh.cc

    
    io/output.cc
    io/output_vtk.cc
    io/output_msh.cc
    io/output_tmp.cc
    
    mesh/msh_gmshreader.cc
    mesh/boundaries.cc
    mesh/edges.cc
    mesh/elements.cc
    mesh/mesh.cc
    mesh/neighbours.cc
    mesh/nodes.cc
    mesh/sides.cc
    mesh/topology.cc

    quadrature/quadrature_lib.cc

    reaction/linear_reaction.cc
    semchem/semchem_interface.cc


    transport/transport.cc
    transport/transport_operator_splitting.cc
    transport/transport_dg.cc
    transport/sources.cc



    equation.cc
    time_marks.cc
    materials.cc
    time_governor.cc
    hc_explicit_sequential.cc
    field_p0.cc
)

set(sources_lib_la
	la/distribution.cc
	la/local_to_global_map.cc
    la/schur.cc
    la/solve.cc
    la/linsys.cc
    la/sparse_graph.cc
)

set(sources_lib_system
    system/math_fce.cc
    system/sys_function_stack.cpp
    system/sys_profiler.cc
    system/system.cc
    system/xio.cc
    system/exceptions.cc
    io/read_ini.cc
)

set(sources_lib_input
    #input/data_tree.cpp
    #//input/Generic_node.cpp
    #//input/Record_node.cpp
    #//input/Value_node.cpp
    #//input/Vector_node.cpp
    
    input/type_base.cc
    input/type_record.cc
    input/storage.cc
    input/json_to_storage.cc
    
    input/input_interface.cc
)

# make separate library for Semchem (God save us!)
set(sources_lib_semchem
    semchem/che_semchem.cc
    semchem/che_read.cc
#    semchem/cti_ichnew.c
)

# this seems as converting *.c to *.cc for Semchem files
# but this can be dengerous
#foreach(file ${libflow123d_C_SRCS})
#  set(libflow123d_SRCS ${libFlow123d_SRCS} ${file}c)
#  list(APPEND libflow123d_SRCS ${file}c)
#  configure_file(${file} ${file}c COPYONLY)
#endforeach(file)

# add build targets
#link_directories(${Flow123dLinkDirectories})


# add individual libraries

# system_lib
add_library(system_lib ${sources_lib_system})
target_link_libraries(system_lib ${PETSC_LIBRARIES} ${Boost_LIBRARY_DIRS})

# input_lib
add_library(input_lib ${sources_lib_input})
target_link_libraries(input_lib 
    system_lib json_spirit)

# la_lib
add_library(la_lib ${sources_lib_la})
target_link_libraries(la_lib 
    input_lib system_lib
    ${PETSC_LIBRARIES} ${Armadillo_LIBRARIES} ${Boost_LIBRARY_DIRS})

# all remaining sources 
add_library(flow123d_lib ${sources_lib_non_system})
target_link_libraries(flow123d_lib 
    input_lib system_lib la_lib fparser semchem
    ${PETSC_LIBRARIES} ${Armadillo_LIBRARIES} ${Boost_LIBRARY_DIRS})
set_target_properties(flow123d_lib PROPERTIES
    ARCHIVE_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/lib
)

# add target for libsemchem
add_library(semchem ${sources_lib_semchem})
set_target_properties(semchem PROPERTIES
  ARCHIVE_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/lib
)







set_directory_properties(PROPERTIES ADDITIONAL_MAKE_CLEAN_FILES "${CMAKE_BINARY_DIR}/bin/current_flow;${CMAKE_BINARY_DIR}/bin/mpiexec")


