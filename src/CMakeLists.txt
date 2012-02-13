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
message(STATUS "Include directories: ${Flow123d_SOURCE_DIR}")

set(libflow123d_SRCS
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
    materials.cc
    msh_gmshreader.cc
    output.cc
    problem.cc
    read_ini.cc    
    solve.cc
    sparse_graph.cc
    time_governor.cc
    transport.cc
    sources.cc
    transport_operator_splitting.cc
    transport_operator_splitting.hh
    xio.cc
)

set(libflow123d_C_SRCS
    semchem/che_semchem.c
    semchem/cti_ichnew.c
)

foreach(file ${libflow123d_C_SRCS})
#  set(libflow123d_SRCS ${libFlow123d_SRCS} ${file}c)
  list(APPEND libflow123d_SRCS ${file}c)
  configure_file(${file} ${file}c COPYONLY)
endforeach(file)

# add build targets
link_directories(${Flow123dLinkDirectories})
add_library(Flow123d ${libflow123d_SRCS})
set_target_properties(Flow123d PROPERTIES
  ARCHIVE_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/lib
)
 
add_executable(flow123d main.cc)
# creates soft link to machine dependent flow runner
if(EXISTS ${CMAKE_SOURCE_DIR}/bin/${MACHINE}_flow.sh)
  add_custom_command(TARGET flow123d POST_BUILD
    COMMAND rm ARGS -f ${CMAKE_BINARY_DIR}/bin/current_flow 
    COMMAND ln ARGS -s ${CMAKE_SOURCE_DIR}/bin/${MACHINE}_flow.sh ${CMAKE_BINARY_DIR}/bin/current_flow
  )
else(EXISTS ${CMAKE_SOURCE_DIR}/bin/${MACHINE}_flow.sh)
  add_custom_command(TARGET flow123d POST_BUILD
    COMMAND rm ARGS -f ${CMAKE_BINARY_DIR}/bin/current_flow 
    COMMAND ln ARGS -s ${CMAKE_SOURCE_DIR}/bin/generic_flow.sh ${CMAKE_BINARY_DIR}/bin/current_flow
  )
endif(EXISTS ${CMAKE_SOURCE_DIR}/bin/${MACHINE}_flow.sh)
set_source_files_properties(${CMAKE_BINARY_DIR}/bin/current_flow PROPERTIES GENERATED true SYMBOLIC true)

#creates soft link to mpiexec
if(EXISTS ${MPIEXEC})
  add_custom_command(TARGET flow123d POST_BUILD
    COMMAND rm ARGS -f ${CMAKE_BINARY_DIR}/bin/mpiexec
    COMMAND ln ARGS -s ${MPIEXEC} ${CMAKE_BINARY_DIR}/bin/mpiexec
  )
  set_source_files_properties(${CMAKE_BINARY_DIR}/bin/mpiexec PROPERTIES GENERATED true SYMBOLIC true)
elseif(EXISTS ${PETSC_DIR}/${PETSC_ARCH}/bin/mpiexec)
  add_custom_command(TARGET flow123d POST_BUILD
    COMMAND rm ARGS -f ${CMAKE_BINARY_DIR}/bin/mpiexec
    COMMAND ln ARGS -s ${PETSC_DIR}/${PETSC_ARCH}/bin/mpiexec ${CMAKE_BINARY_DIR}/bin/mpiexec
  )
  set_source_files_properties(${CMAKE_BINARY_DIR}/bin/mpiexec PROPERTIES GENERATED true SYMBOLIC true)
else()
  message(WARNING "Softlink to mpiexec wasn't created.")
endif()
set_directory_properties(PROPERTIES ADDITIONAL_MAKE_CLEAN_FILES "${CMAKE_BINARY_DIR}/bin/current_flow;${CMAKE_BINARY_DIR}/bin/mpiexec")


target_link_libraries(flow123d ${PETSC_LIBRARIES} ${PARMETIS_LIBRARIES} Flow123d)
set_target_properties(flow123d PROPERTIES
   RUNTIME_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/bin
   OBJECT_DIR ${PROJECT_BINARY_DIR}/build/obj
)