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
#set(USE_DEAL_II true)
if (USE_DEAL_II)
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
    unset(DEAL_II_COMPONENT_PACKAGE  CACHE)

    set(DEAL_II_COMPONENT_EXAMPLES OFF)
    set(DEAL_II_COMPONENT_COMPAT_FILES OFF)
    set(DEAL_II_COMPONENT_DOCUMENTATION OFF)
    set(DEAL_II_COMPONENT_MESH_CONVERTER OFF)
    set(DEAL_II_COMPONENT_PARAMETER_GUI OFF)

#add_subdirectory("dealii")
    set(DEAL_II_LIBRARY_NAME ${DEAL_II_LIBRARY_NAME} PARENT_SCOPE)
    message(STATUS "===========================================")


endif()


# patched json_spirit lib
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
    system/exceptions.cc
    system/stack_trace.cc
    system/asserts.cc
    system/file_path.cc
    system/tokenizer.cc
    system/application_base.cc
    system/logger.cc
    system/logger_options.cc
    system/armadillo_tools.cc
    system/fmt/format.cc
    system/fmt/ostream.cc
    system/fmt/posix.cc
)
target_link_libraries(system_lib 
    ${PETSC_LIBRARIES}  
    ${PETSC_EXTERNAL_LIB} 
    ${Boost_LIBRARIES} 
    ${PYTHON_LIBRARIES}
    nlohmann_json::nlohmann_json)
add_dependencies(system_lib update_source_tree)


# input_lib
add_library(input_lib 
    input/type_base.cc
    input/type_record.cc
    input/type_abstract.cc
    input/type_selection.cc
    input/type_generic.cc
    input/type_tuple.cc
    input/storage.cc
    input/reader_to_storage.cc
    input/reader_internal_base.cc
    input/reader_internal.cc
    input/reader_internal_transpose.cc
    input/reader_internal_csv.cc
    input/path_base.cc
    input/path_json.cc
    input/path_yaml.cc
    input/type_output.cc
    
    input/accessors.cc
    
    input/csv_tokenizer.cc
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
    la/bddcml_wrapper.cc
    la/linsys_BDDC.cc
    la/linsys_PETSC.cc
    la/sparse_graph.cc
    la/local_system.cc
    la/vector_mpi.cc
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
    tools/unit_si.cc
    tools/unit_converter.cc
)
target_link_libraries(tools_lib
    system_lib input_lib )



# mesh_lib
add_library(mesh_lib 
    mesh/region.cc
    mesh/elements.cc
    mesh/duplicate_nodes.cc
    mesh/mesh.cc
    mesh/bc_mesh.cc
    mesh/neighbours.cc
    mesh/accessors.cc
#    mesh/intersection.cc
    mesh/ref_element.cc
    mesh/partitioning.cc
    mesh/region_set.cc
    
    mesh/bounding_box.cc
    mesh/bih_tree.cc
    
#     mesh/ngh/src/abscissa.cpp
#     mesh/ngh/src/bisector.cpp
#     mesh/ngh/src/intersection.cpp
#     mesh/ngh/src/intersectionLocal.cpp
#     mesh/ngh/src/mathfce.cpp
#     mesh/ngh/src/matrix.cpp
#     mesh/ngh/src/plain.cpp
#     mesh/ngh/src/point.cpp
#     mesh/ngh/src/polygon.cpp
#     mesh/ngh/src/tetrahedron.cpp
#     mesh/ngh/src/triangle.cpp
#     mesh/ngh/src/vector.cpp
#     mesh/ngh/src/vertex.cpp

    intersection/plucker.cc
    intersection/intersection_point_aux.cc
    intersection/intersection_aux.cc
    intersection/intersection_local.cc
    intersection/compute_intersection.cc
    intersection/inspect_elements_algorithm.cc
    intersection/mixed_mesh_intersections.cc
)

target_link_libraries(mesh_lib
    input_lib system_lib la_lib
    armadillo 
    ${Boost_LIBRARIES}
    ${PugiXml_LIBRARY}
    ${Zlib_LIBRARY})

# fem_lib    
add_library(fem_lib 
    quadrature/quadrature.cc
    quadrature/quadrature_lib.cc
    quadrature/intersection_quadrature.cc
    fem/discrete_space.cc
    fem/dofhandler.cc
    fem/finite_element.cc
    fem/fe_p.cc
    fem/fe_rt.cc
    fem/fe_system.cc
    fem/element_values.cc
    fem/fe_values.cc
    fem/fe_values_views.cc
    fem/mapping_p1.cc
)
target_link_libraries(fem_lib
    system_lib mesh_lib 
    armadillo 
    ${Boost_LIBRARIES})



# io_lib  
add_library(io_lib 
    io/msh_basereader.cc
    io/msh_gmshreader.cc
    io/msh_vtkreader.cc
    io/msh_pvdreader.cc
    io/element_data_cache.cc
    io/reader_cache.cc

    io/output_time.cc
    io/output_vtk.cc
    io/output_msh.cc
    io/observe.cc
    io/output_mesh.cc
    io/output_time_set.cc
)

target_link_libraries(io_lib
    system_lib input_lib mesh_lib fem_lib la_lib
    armadillo 
    ${Boost_LIBRARIES}
    ${PugiXml_LIBRARY}
    ${Zlib_LIBRARY})





# coupling_lib  
add_library(coupling_lib 
   # unfortunatelly field sources depends on output sources
   # and viceversa
    
    fields/equation_output.cc
    fields/field_values.cc
    fields/field_common.cc
    fields/field.cc
    fields/bc_field.cc
    fields/multi_field.cc
    fields/bc_multi_field.cc
    fields/field_algo_base.cc    
    fields/field_flag.cc
    fields/field_set.cc
    fields/generic_field.cc
    fields/field_constant.cc
    fields/field_formula.cc
    # fields/field_elementwise.cc
    # fields/field_interpolated_p0.cc
    fields/table_function.cc
    fields/field_time_function.cc
    fields/field_divide.cc
    fields/field_fe.cc
    fields/fe_value_handler.cc
    fields/eval_subset.cc
    fields/eval_points.cc
    fields/field_value_cache.cc
    fields/surface_depth.cc
    coupling/equation.cc
    coupling/balance.cc
    # coupling/hc_explicit_sequential.cc
)
target_link_libraries(coupling_lib
    system_lib input_lib mesh_lib tools_lib fem_lib io_lib
    fparser     
    armadillo
    ${Boost_LIBRARIES})



# make separate library for Semchem (God save us!)
#add_library(semchem 
#    semchem/che_semchem.cc
#    semchem/che_read.cc
#)
#set_target_properties(semchem 
#    PROPERTIES ARCHIVE_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/lib
#)



set(sources_lib_non_system
)
# all remaining sources 
add_library(flow123d_lib 
    flow/darcy_flow_mh_output.cc
    flow/darcy_flow_mh.cc
    flow/darcy_flow_lmh.cc
    flow/soil_models.cc
    flow/richards_lmh.cc
    flow/mortar_assembly.cc

    reaction/reaction_term.cc
    reaction/first_order_reaction.cc
    reaction/first_order_reaction_base.cc
    reaction/radioactive_decay.cc
    reaction/sorption_base.cc
    reaction/sorption.cc
    reaction/dual_porosity.cc
    reaction/isotherm.cc
    reaction/linear_ode_solver.cc
    #semchem/semchem_interface.cc

    transport/concentration_model.cc
    transport/heat_model.cc
    transport/transport.cc
    transport/transport_operator_splitting.cc
    transport/transport_dg.cc
    
    transport/substance.cc

    mechanics/elasticity.cc

    coupling/hc_explicit_sequential.cc
    coupling/hm_iterative.cc
)    
    
target_link_libraries(flow123d_lib 
    fem_lib mesh_lib la_lib input_lib io_lib
    system_lib tools_lib coupling_lib
    fparser  
    armadillo 
#    ${TBB_LIBRARIES}
#    ${DEAL_II_LIBRARY_NAME}
    ${PETSC_LIBRARIES} 
    ${BDDCML_LIBRARIES} 
    ${PETSC_EXTERNAL_LIB} 
    ${Boost_LIBRARIES}
)
set_target_properties(flow123d_lib 
    PROPERTIES ARCHIVE_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/lib
)
#add_dependencies(flow123d_lib tbb)

