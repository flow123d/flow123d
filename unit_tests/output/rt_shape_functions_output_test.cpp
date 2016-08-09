/*
 * rt_shape_functions_output_test.cpp
 *
 *  Created on: July 1, 2016
 *      Author: pe
 */

#define TEST_USE_PETSC
#define FEAL_OVERRIDE_ASSERTS
#include <flow_gtest_mpi.hh>
#include "io/output_time.hh"
#include "io/output_vtk.hh"
#include "mesh/mesh.h"
#include "input/reader_to_storage.hh"
#include "system/sys_profiler.hh"

#include <fields/field_fe.hh>
#include <fields/field_set.hh>
#include <fields/field_common.hh>
#include "input/input_type.hh"

#include "io/output_mesh.hh"

FLOW123D_FORCE_LINK_IN_PARENT(field_constant)

const string input_rt = R"INPUT(
{   
   output_stream = {
    file = "./test_shape.pvd", 
    format = {
        TYPE = "vtk", 
        variant = "ascii"
    },
    output_mesh = {
        max_level = 6
    }
  }
}
)INPUT";

// simplest mesh
string ref_element_mesh = R"CODE(
$MeshFormat
2.2 0 8
$EndMeshFormat
$Nodes
3
1 0 0 0
2 1 0 0
3 0 1 0
$EndNodes
$Elements
1
1 2 2 39 40 1 2 3
$EndElements
)CODE";

#include "fem/mapping_p1.hh"
#include "fem/fe_rt.hh"
#include "fem/dofhandler.hh"
#include "fem/fe_values.hh"
#include "quadrature/quadrature_lib.hh"
#include "fields/field_fe.hh"

#include "fem/singularity.hh"
#include "fem/fe_rt_xfem.hh"
#include "fem/fe_p0_xfem.hh"

void output_field_fe(FiniteElement<1,3>& fe_1,
                     FiniteElement<2,3>& fe_2,
                     FiniteElement<3,3>& fe_3,
                     const std::map<unsigned int, double>& dof_values,
                     bool is_scalar)
{
// read mesh - simplset cube from test1
    Mesh* mesh = new Mesh();
    stringstream in(ref_element_mesh.c_str());
    mesh->read_gmsh_from_stream(in);
    
    DOFHandlerMultiDim dofhandler(*mesh);
       
    MappingP1<1,3> map1;
    MappingP1<2,3> map2;
    MappingP1<3,3> map3;
    
    dofhandler.distribute_dofs(fe_1, fe_2, fe_3);
    
    Vec data_vec;
    VecCreateSeq(PETSC_COMM_SELF, dofhandler.n_global_dofs(), &data_vec);
    
    for(auto &pair: dof_values){
        VecSetValue(data_vec, pair.first, pair.second, ADD_VALUES);
    }
    
    FieldCommon* output_field;
        
    if(is_scalar) {
        std::shared_ptr<FieldFE<3,FieldValue<3>::Scalar>> field_fe = std::make_shared<FieldFE<3,FieldValue<3>::Scalar>>(1);
        field_fe->set_mesh(mesh,false);
        field_fe->set_fe_data(&dofhandler,&map1, &map2, &map3, &data_vec);
    
        Field<3,FieldValue<3>::Scalar>* of = new Field<3,FieldValue<3>::Scalar>();
        of->set_mesh(*mesh);
        of->set_field(mesh->region_db().get_region_set("ALL"), field_fe, 0);
        output_field = of;
    }
    else {
        std::shared_ptr<FieldFE<3,FieldValue<3>::VectorFixed>> field_fe 
                        = std::make_shared<FieldFE<3,FieldValue<3>::VectorFixed>>(1);
        field_fe->set_mesh(mesh,false);
        field_fe->set_fe_data(&dofhandler,&map1, &map2, &map3, &data_vec);
        
        Field<3,FieldValue<3>::VectorFixed>* of = new Field<3,FieldValue<3>::VectorFixed>();
        of->set_mesh(*mesh);
        of->set_field(mesh->region_db().get_region_set("ALL"), field_fe, 0);
        output_field = of;
    }
    output_field->output_type(OutputTime::CORNER_DATA);
    
    // create field set of output fields
    FieldSet output_fields;
    output_fields += output_field->name("shape_func").units(UnitSI::dimensionless()).flags_add(FieldFlag::allow_output);
    
    // set time to all fields
    output_fields.set_time(0.0, LimitSide::right);

//     FilePath::set_io_dirs(".",UNIT_TESTS_SRC_DIR,"",".");
    FilePath::set_io_dirs(".",UNIT_TESTS_SRC_DIR,"","output");
    // read input string
    Input::Record in_rec = Input::ReaderToStorage( input_rt,
                                   Input::Type::Record("ErrorFieldTest","")
                                   .declare_key("output_stream", OutputTime::get_input_type(), Input::Type::Default::obligatory(), "")
                                   .close(),
                                   Input::FileFormat::format_JSON ).get_root_interface<Input::Record>();
                                   
    // create output
    std::shared_ptr<OutputTime> output = std::make_shared<OutputVTK>();
    output->init_from_input("dummy_equation", in_rec.val<Input::Record>("output_stream"));
    
    // register output fields, compute and write data
    output_fields.output(output);
    output->write_time_frame();
    
    delete output_field;
    delete mesh;
}


TEST(ShapeFunctionOutput, rt_shape) {
    
    FE_RT0<1,3> fe_rt1;
    FE_RT0<2,3> fe_rt2;
    FE_RT0<3,3> fe_rt3;
    
    std::map<unsigned int, double> dof_values = {
        { 0, 1.0 }
//         { 1, 1.0 }
//         { 2, 1.0 }
    };
    
    output_field_fe(fe_rt1, fe_rt2, fe_rt3, dof_values, false);
}


TEST(ShapeFunctionOutput, rt_xfem_shape) {
    
    Singularity0D<3> func({0.2,0.2,0},0.05);
    
    FE_RT0<1,3> fe_rt1;
    FE_RT0<2,3> fe_rt2;
    FE_RT0_XFEM<2,3> fe_rt_xfem(&fe_rt2,{&func});
    FE_RT0<3,3> fe_rt3;
    
//     //precise enrichment function approx.
//     std::map<unsigned int, double> dof_values = {
//         { 0, 1.53846153846154 },    // interpolation dofs
//         { 1, 1.53846153846154 },
//         { 2, 2.35702260395516 },
//         { 3, 1.0 },
//         { 4, 1.0 },
//         { 5, 1.0 }
//     };
    
    std::map<unsigned int, double> dof_values = {
//        { 3, 1.0 }
//        { 4, 1.0 }
       { 5, 1.0 }
    };
    
    output_field_fe(fe_rt1, fe_rt_xfem, fe_rt3, dof_values, false);
}


TEST(ShapeFunctionOutput, p0_xfem) {

    Singularity0D<3> func({0.2,0.2,0},0.05);
    
    FE_P_disc<0,1,3> fe_p_1;
    FE_P_disc<0,2,3> fe_p_2;
    FE_P0_XFEM<2,3> fe_p0_xfem(&fe_p_2,{&func});
    FE_P_disc<0,3,3> fe_p_3;

//     //precise enrichment function approx.    
    std::map<unsigned int, double> dof_values = {
       { 0, -0.192831240405992 },   // interpolation dof
       { 1, 1.0 },
       { 2, 1.0 },
       { 3, 1.0 }
    };
    
//         std::map<unsigned int, double> dof_values = {
// //        { 0, 1.0 }
// //        { 1, 1.0 }
// //        { 2, 1.0 }
//        { 3, 1.0 }
//     };
    
    output_field_fe(fe_p_1, fe_p0_xfem, fe_p_3, dof_values, true);
}