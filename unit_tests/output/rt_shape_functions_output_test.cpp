/*
 * rt_shape_functions_output_test.cpp
 *
 *  Created on: July 1, 2016
 *      Author: pe
 */

#define TEST_USE_PETSC
#define FEAL_OVERRIDE_ASSERTS
#include <flow_gtest_mpi.hh>

#include "fields/equation_output.hh"
#include "io/output_time.hh"
#include "io/output_vtk.hh"
#include "io/output_mesh.hh"

#include "mesh/mesh.h"
#include "mesh_constructor.hh"
#include "input/reader_to_storage.hh"
#include "input/input_type.hh"
#include "system/sys_profiler.hh"

#include <fields/field_fe.hh>
#include <fields/field_set.hh>
#include <fields/field_common.hh>

#include "fem/mapping_p1.hh"
#include "quadrature/quadrature_lib.hh"
#include "fem/dofhandler.hh"
#include "fem/fe_values.hh"

#include "fem/fe_rt.hh"
#include "fem/singularity.hh"
#include "fem/fe_rt_xfem.hh"
#include "fem/fe_rt_xfem_single.hh"
#include "fem/fe_p0_xfem.hh"



FLOW123D_FORCE_LINK_IN_PARENT(field_constant)

const std::string input_rt = R"INPUT(
{   
   output_stream = {
    file = "test_shape", 
    format = {
        TYPE = "vtk", 
        variant = "ascii"
    },
    output_mesh = {
        max_level = 5
    }
  }
  ,output = {fields = ["shape_func"]}
}
)INPUT";

const std::string mesh_file_2d = "ref_element_mesh_2d.msh";
const std::string mesh_file_3d = "ref_element_mesh_3d.msh";

bool replace_string(std::string& str, const std::string& from, const std::string& to) {
    size_t start_pos = str.find(from);
    if(start_pos == std::string::npos)
        return false;
    str.replace(start_pos, from.length(), to);
    return true;
}


void output_field_fe(FiniteElement<1,3>& fe_1,
                     FiniteElement<2,3>& fe_2,
                     FiniteElement<3,3>& fe_3,
                     const std::map<unsigned int, double>& dof_values,
                     bool is_scalar,
                     std::string file_name,
                     std::string input_mesh_file)
{
    // replace correct output file name in input record string
    std::string input_json = input_rt;
    bool res = replace_string(input_json, "test_shape" ,file_name);
    ASSERT(res);
    
    // setup FilePath directories
    FilePath::set_io_dirs(".",FilePath::get_absolute_working_dir(),"",".");
    // read mesh
    Mesh *mesh = mesh_full_constructor("{mesh_file=\""+input_mesh_file+"\"}");
    
    shared_ptr<DOFHandlerMultiDim> dofhandler = std::make_shared<DOFHandlerMultiDim>(*mesh);
       
    MappingP1<1,3> map1;
    MappingP1<2,3> map2;
    MappingP1<3,3> map3;
    
    dofhandler->distribute_dofs(fe_1, fe_2, fe_3);
    
    VectorSeqDouble vvv;
    vvv.resize(dofhandler->n_global_dofs());
    for(auto &pair: dof_values){
        VecSetValue(vvv.get_data_petsc(), pair.first, pair.second, ADD_VALUES);
    }
    
    FieldCommon* output_field;
    
    if(is_scalar) {
        std::shared_ptr<FieldFE<3,FieldValue<3>::Scalar>> field_fe = std::make_shared<FieldFE<3,FieldValue<3>::Scalar>>(1);
        field_fe->set_mesh(mesh,false);
        field_fe->set_fe_data(dofhandler,&map1, &map2, &map3, &vvv);
    
        Field<3,FieldValue<3>::Scalar>* of = new Field<3,FieldValue<3>::Scalar>();
        of->set_mesh(*mesh);
        of->set_field(mesh->region_db().get_region_set("ALL"), field_fe, 0);
        output_field = of;
    }
    else {
        std::shared_ptr<FieldFE<3,FieldValue<3>::VectorFixed>> field_fe 
                        = std::make_shared<FieldFE<3,FieldValue<3>::VectorFixed>>(1);
        field_fe->set_mesh(mesh,false);
        field_fe->set_fe_data(dofhandler,&map1, &map2, &map3, &vvv);
        
        Field<3,FieldValue<3>::VectorFixed>* of = new Field<3,FieldValue<3>::VectorFixed>();
        of->set_mesh(*mesh);
        of->set_field(mesh->region_db().get_region_set("ALL"), field_fe, 0);
        output_field = of;
    }
    output_field->output_type(OutputTime::CORNER_DATA);
    
    // create field set of output fields
    EquationOutput output_fields;
    output_fields += output_field->name("shape_func").units(UnitSI::dimensionless()).flags_add(FieldFlag::allow_output);
    
    // set time to all fields
    output_fields.set_time(0.0, LimitSide::right);

//     FilePath::set_io_dirs(".",UNIT_TESTS_SRC_DIR,"",".");
    FilePath::set_io_dirs(".",UNIT_TESTS_SRC_DIR,"","output");
    
    // create input record
    Input::Type::Record rec_type = Input::Type::Record("ErrorFieldTest","")
        .declare_key("output_stream", OutputTime::get_input_type(), Input::Type::Default::obligatory(), "")
        .declare_key("output", output_fields.make_output_type("test_eq"), Input::Type::Default::obligatory(), "")
        .close();

    // read input string
    Input::ReaderToStorage reader( input_json, rec_type, Input::FileFormat::format_JSON );
    Input::Record in_rec=reader.get_root_interface<Input::Record>();
    
        
    // create output
    std::shared_ptr<OutputTime> output = std::make_shared<OutputVTK>();
    output->init_from_input("dummy_equation", *mesh, in_rec.val<Input::Record>("output_stream"));
    output_fields.initialize(output, in_rec.val<Input::Record>("output"), TimeGovernor());
    
    // register output fields, compute and write data
    output_fields.output(0.0);
    output->write_time_frame();
    
    delete output_field;
    delete mesh;
}


TEST(ShapeFunctionOutput, rt_xfem_shape_2d) {
    
    auto func = std::make_shared<Singularity<0>>(arma::vec({0.2,0.2,0}),0.05,
                                                 arma::vec({0, 0, 1}), arma::vec({0, 0, 1}), 100);
    
    FE_RT0<1,3> fe_rt1;
    FE_RT0<2,3> fe_rt2;
    FE_RT0_XFEM_S<2,3> fe_rt_xfem(&fe_rt2,{func});
    FE_RT0<3,3> fe_rt3;
    
   std::map<unsigned int, double> dof_values;

    // print all shape functions
    std::string filename = "test_rt_2d_";
    for(unsigned int i=0; i < fe_rt_xfem.n_dofs(); i++){
        dof_values = {{i, 1.0}};  //select only one shape function
        output_field_fe(fe_rt1, fe_rt_xfem, fe_rt3, dof_values, false, filename + std::to_string(i), mesh_file_2d);
    }
    
//     //     //precise enrichment function approx.
//     dof_values = {
//         { 0, 1.53846153846154 },    // interpolation dofs
//         { 1, 1.53846153846154 },
//         { 2, 2.35702260395516 },
//         { 3, 1.0 },
//         { 4, 1.0 },
//         { 5, 1.0 }
//     };

    //     //precise enrichment function approx.
//     dof_values = {
//         { 0, 2.111270793766706 },    // interpolation dofs
//         { 1, 2.111270793766706 },
//         { 2, 2.060775278298543 },
//         { 3, 1.0 }
        dof_values = {
        { 0, -0.3360070022065689 },    // interpolation dofs
        { 1, -0.33600700220656896 },
        { 2, -0.3279799752702353 },
        { 3, 1.0 }
    };
    
    output_field_fe(fe_rt1, fe_rt_xfem, fe_rt3, dof_values, false, "test_rt_2d",  mesh_file_2d);
}


TEST(ShapeFunctionOutput, rt_xfem_shape_3d) {
    
    auto func = std::make_shared<Singularity<1>>(arma::vec({0.2,0.2,0}),arma::vec({0.2,0.2,1.0}),0.05, 50, 20);
    
    FE_RT0<1,3> fe_rt1;
    FE_RT0<2,3> fe_rt2;
    FE_RT0<3,3> fe_rt3;
    FE_RT0_XFEM_S<3,3> fe_rt_xfem(&fe_rt3,{func});
    
   std::map<unsigned int, double> dof_values;

    // print all shape functions
    std::string filename = "test_rt__3d_";
//     for(unsigned int i=0; i < fe_rt_xfem.n_dofs(); i++){
//         dof_values = {{i, 1.0}};  //select only one shape function
//         output_field_fe(fe_rt1, fe_rt2, fe_rt_xfem, dof_values, false, filename + std::to_string(i), mesh_file_3d);
//     }
    
//     //     //precise enrichment function approx.
//     dof_values = {
//         { 0, 1.53846153846154 },    // interpolation dofs
//         { 1, 1.53846153846154 },
//         { 2, 2.35702260395516 },
//         { 3, 1.0 },
//         { 4, 1.0 },
//         { 5, 1.0 }
//     };

    //     //precise enrichment function approx.
//     dof_values = {
//         { 0, 2.111270793766706 },    // interpolation dofs
//         { 1, 2.111270793766706 },
//         { 2, 2.060775278298543 },
//         { 3, 1.0 }
        dof_values = {
        { 0, 0.0 },    // interpolation dofs
        { 1, -0.1173743870500195 },    // interpolation dofs
        { 2, -0.1173743870500195 },
        { 3, -0.06525196651968521 },
        { 4, 1.0 }
    };
    
    output_field_fe(fe_rt1, fe_rt2, fe_rt_xfem, dof_values, false, "test_rt_3d", mesh_file_3d);
}

// TEST(ShapeFunctionOutput, p0_xfem) {
// 
//     auto func = std::make_shared<Singularity0D>(arma::vec({0.2,0.2,0}),0.05,
//                                                    arma::vec({0, 0, 1}), arma::vec({0, 0, 1}));
//     
//     FE_P_disc<0,1,3> fe_p_1;
//     FE_P_disc<0,2,3> fe_p_2;
//     FE_P0_XFEM<2,3> fe_p0_xfem(&fe_p_2,{func});
//     FE_P_disc<0,3,3> fe_p_3;
// 
//     std::map<unsigned int, double> dof_values;
// 
//     // print all shape functions
//     std::string filename = "test_p0_";
//     for(unsigned int i=0; i < fe_p0_xfem.n_dofs(); i++){
//         dof_values = {{i, 1.0}};  //select only one shape function
//         output_field_fe(fe_p_1, fe_p0_xfem, fe_p_3, dof_values, true, filename + std::to_string(i));
//     }
//     
// //     //precise enrichment function approx.    
//     dof_values = {
//        { 0, -0.192831240405992 },   // value of enrich function at interpolation dof point
//        { 1, 1.0 },
//        { 2, 1.0 },
//        { 3, 1.0 }
//     };
//     
//     output_field_fe(fe_p_1, fe_p0_xfem, fe_p_3, dof_values, true, "test_p0");
// }