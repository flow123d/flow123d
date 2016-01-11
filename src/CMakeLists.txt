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


# add deal.II library
message(STATUS "===========================================")
# Note: deal.II implicitly uses variables PETSC_DIR, PETSC_ARCH from upstream
set(DEAL_II_WITH_MUPARSER OFF)
set(DEAL_II_WITH_MPI ON)
set(DEAL_II_WITH_THREADS OFF)
set(BOOST_DIR "${BOOST_ROOT}")
set(BOOST_LIBRARIES "${BOOST_LIBRARY_DIRS}")
#set(TBB_FOUND "true")
set(DEAL_II_WITH_UMFPACK OFF)
## Check if UMFPACK is installed with PETSc
#if (PETSC_VAR_UMFPACK_LIB)
#  set(UMFPACK_FOUND "true")
#  set(UMFPACK_LINKER_FLAGS ${PETSC_VAR_UMFPACK_LIB})
#else()
#  MESSAGE(FATAL_ERROR "Deal.II depends on UMFPACK, but PETSc is not built with it.")
##  set(DEAL_II_ALLOW_BUNDLED OFF)
##  set(DEAL_II_FORCE_BUNDLED_UMFPACK ON)
#endif()
set(DEAL_II_WITH_LAPACK ON)
set(LAPACK_FOUND "true")
set(LAPACK_LIBRARIES "${PETSC_EXTERNAL_LIB}")
set(LAPACK_LINKER_FLAGS ${PETSC_VAR_BLASLAPACK_LIB})
add_subdirectory("dealii")
set(DEAL_II_LIBRARY_NAME ${DEAL_II_LIBRARY_NAME} PARENT_SCOPE)
message(STATUS "===========================================")


add_subdirectory("input/json_spirit")

set(CMAKE_INCLUDE_CURRENT_DIR ON)
#message(STATUS "Include directories: ${flow123d_includes}")
#include_directories(${flow123d_includes})


# add individual libraries and object groups

if (mpiuni_sources) 
  add_library(mpiuni_lib ${mpiuni_sources})
endif()

# when the precompiled headers are used in the future
#if (${PRECOMPILED_HEADERS})
#    ADD_PRECOMPILED_HEADER(system_lib "pre_system.hh")
#endif()

# system_lib
add_library(system_lib 
    system/python_loader.cc
    system/math_fce.cc
    system/sys_profiler.cc
    system/time_point.cc
    system/system.cc
    system/xio.cc
    system/exceptions.cc
    io/read_ini.cc
    system/file_path.cc
    system/tokenizer.cc
    system/application_base.cc
)
target_link_libraries(system_lib 
    ${PETSC_LIBRARIES}  
    ${PETSC_EXTERNAL_LIB} 
    ${Boost_LIBRARIES} 
    ${PYTHON_LIBRARIES})
add_dependencies(system_lib update_source_tree)


# input_lib
add_library(input_lib 
    input/type_base.cc
    input/type_record.cc
    input/type_abstract.cc
    input/type_selection.cc
    input/type_generic.cc
    input/storage.cc
    input/reader_to_storage.cc
    input/path_base.cc
    input/path_json.cc
    input/path_yaml.cc
    input/type_output.cc
    
    input/accessors.cc
    input/storage_transpose.cc
)
target_link_libraries(input_lib 
    system_lib json_spirit
    ${YamlCpp_LIBRARY})


# la_lib
add_library(la_lib 
    la/distribution.cc
    la/local_to_global_map.cc
    la/schur.cc
    la/linsys.cc
    la/linsys_BDDC.cc
    la/linsys_PETSC.cc
    la/sparse_graph.cc
)
target_link_libraries(la_lib 
    input_lib system_lib
    armadillo 
    ${PETSC_LIBRARIES}  
    ${BDDCML_LIBRARIES} 
    ${PETSC_EXTERNAL_LIB}  
    ${Boost_LIBRARIES})
set_target_properties(la_lib 
    PROPERTIES  COMPILE_FLAGS "${BDDCML_CDEFS}" ) 



# tools_lib    
add_library(tools_lib 
    tools/interpolant.cc
    tools/adaptivesimpson.cc
    tools/time_marks.cc
    tools/time_governor.cc
)
target_link_libraries(tools_lib
    system_lib input_lib )



# mesh_lib
add_library(mesh_lib 
        mesh/region.cc
    mesh/msh_gmshreader.cc
    mesh/boundaries.cc
    mesh/edges.cc
    mesh/elements.cc
    mesh/mesh.cc
    mesh/neighbours.cc
    mesh/sides.cc
    mesh/intersection.cc
    mesh/ref_element.cc
    mesh/partitioning.cc
    
    mesh/bounding_box.cc
    mesh/bih_tree.cc
    mesh/reader_instances.cc
    
    mesh/ngh/src/abscissa.cpp
    mesh/ngh/src/bisector.cpp
    mesh/ngh/src/intersection.cpp
    mesh/ngh/src/intersectionLocal.cpp
    mesh/ngh/src/mathfce.cpp
    mesh/ngh/src/matrix.cpp
    mesh/ngh/src/plain.cpp
    mesh/ngh/src/point.cpp
    mesh/ngh/src/polygon.cpp
    mesh/ngh/src/tetrahedron.cpp
    mesh/ngh/src/triangle.cpp
    mesh/ngh/src/vector.cpp
    mesh/ngh/src/vertex.cpp
)
target_link_libraries(mesh_lib
    input_lib system_lib la_lib 
    armadillo 
    ${Boost_LIBRARIES})



# fem_lib    
add_library(fem_lib 
    quadrature/quadrature_lib.cc
    fem/fe_p.cc
    fem/dofhandler.cc
    fem/fe_values.cc
    fem/finite_element.cc
    fem/mapping_p1.cc
)
target_link_libraries(fem_lib
    system_lib mesh_lib 
    armadillo 
    ${Boost_LIBRARIES})






# coupling_lib  
add_library(coupling_lib 
   # unfortunatelly field sources depends on output sources
   # and viceversa
    io/output_time.cc
    io/output_vtk.cc
    io/output_msh.cc
    
    fields/field_values.cc
    fields/field_common.cc
    fields/field.cc
    fields/bc_field.cc
    fields/multi_field.cc
    fields/field_algo_base.cc    
    fields/field_flag.cc
    fields/field_set.cc
    fields/unit_si.cc
    fields/generic_field.cc
    coupling/equation.cc
    coupling/balance.cc
    # coupling/hc_explicit_sequential.cc
)
target_link_libraries(coupling_lib
    system_lib input_lib mesh_lib tools_lib
    fparser     
    armadillo
    ${Boost_LIBRARIES})



# make separate library for Semchem (God save us!)
add_library(semchem 
    semchem/che_semchem.cc
    semchem/che_read.cc
)
set_target_properties(semchem 
    PROPERTIES ARCHIVE_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/lib
)



set(sources_lib_non_system
)
# all remaining sources 
add_library(flow123d_lib 
    flow/mh_dofhandler.cc
    flow/darcy_flow_mh_output.cc
    flow/darcy_flow_mh.cc
    flow/richards_lmh.cc

    reaction/reaction_term.cc
    reaction/first_order_reaction.cc
    reaction/first_order_reaction_base.cc
    reaction/radioactive_decay.cc
    reaction/pade_approximant.cc
    reaction/sorption_base.cc
    reaction/sorption.cc
    reaction/dual_porosity.cc
    reaction/isotherm.cc
    reaction/linear_ode_solver.cc
    reaction/linear_ode_analytic.cc
    semchem/semchem_interface.cc

    transport/concentration_model.cc
    transport/heat_model.cc
    transport/transport.cc
    transport/transport_operator_splitting.cc
    transport/transport_dg.cc
    
    transport/substance.cc

    coupling/hc_explicit_sequential.cc
)    
target_link_libraries(flow123d_lib 
    fem_lib mesh_lib la_lib semchem  input_lib 
    system_lib tools_lib coupling_lib  tools_lib
    fparser  
    armadillo 
#    ${TBB_LIBRARIES}
    ${DEAL_II_LIBRARY_NAME}
    ${PETSC_LIBRARIES} 
    ${BDDCML_LIBRARIES} 
    ${PETSC_EXTERNAL_LIB} 
    ${Boost_LIBRARIES}
)
set_target_properties(flow123d_lib 
    PROPERTIES ARCHIVE_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/lib
)
add_dependencies(flow123d_lib tbb)

