/*
 * output_vtk_test.cpp
 *
 *  Created on: June 9, 2016
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

#include "fields/field_constant.hh"
#include <fields/field.hh>
#include <fields/field_set.hh>
#include <io/equation_output.hh>
#include <fields/field_common.hh>
#include "input/input_type.hh"

#include "input/accessors.hh"

#include "io/output_mesh.hh"
#include "io/output_element.hh"

FLOW123D_FORCE_LINK_IN_PARENT(field_constant)



TEST(OutputMesh, create_identical)
{
    // setup FilePath directories
    FilePath::set_io_dirs(".",UNIT_TESTS_SRC_DIR,"",".");

    // read mesh - simplset cube from test1
    FilePath mesh_file( string(UNIT_TESTS_SRC_DIR) + "/mesh/simplest_cube.msh", FilePath::input_file);
    Mesh* mesh = new Mesh();
    ifstream in(string(mesh_file).c_str());
    mesh->read_gmsh_from_stream(in);
    
    auto output_mesh = std::make_shared<OutputMesh>(*mesh);
    output_mesh->create_identical_mesh();
    
    std::cout << "nodes: ";
    output_mesh->nodes_->print_all(std::cout);
    std::cout << endl;
    std::cout << "connectivity: ";
    output_mesh->connectivity_->print_all(std::cout);
    std::cout << endl;
    std::cout << "offsets: ";
    output_mesh->offsets_->print_all(std::cout);
    std::cout << endl;
    
    EXPECT_EQ(output_mesh->nodes_->n_values, mesh->n_nodes());
    EXPECT_EQ(output_mesh->offsets_->n_values, mesh->n_elements());
    
    for(const auto &ele : *output_mesh)
    {
    	std::cout << ele.idx() << " " << ele.dim() << "D n_" << ele.n_nodes() << " |";
        for(unsigned int i=0; i < ele.n_nodes(); i++)
        {
        	std::cout << " " << ele.node_index(i);
        }
        std::cout << " |";
        for(auto& v : ele.vertex_list())
        {
        	std::cout << " " << v[0] << " " << v[1] << " " << v[2] << " #";
        }
        std::cout << endl;
        
        ElementAccessor<3> ele_acc = ele.element_accessor();
        EXPECT_EQ(ele.dim(), ele_acc.dim());
        EXPECT_EQ(ele.centre()[0], ele_acc.centre()[0]);
        EXPECT_EQ(ele.centre()[1], ele_acc.centre()[1]);
        EXPECT_EQ(ele.centre()[2], ele_acc.centre()[2]);
    }
    
    auto output_mesh_discont = std::make_shared<OutputMeshDiscontinuous>(*mesh);
    output_mesh_discont->create_mesh(output_mesh);
    
    MessageOut() << "DISCONTINUOUS\n";
    for(const auto &ele : *output_mesh_discont)
    {
    	std::cout << ele.idx() << " " << ele.dim() << "D n_" << ele.n_nodes() << " |";
        for(unsigned int i=0; i < ele.n_nodes(); i++)
        {
        	std::cout << " " << ele.node_index(i);
        }
        std::cout << " |";
        for(auto& v : ele.vertex_list())
        {
        	std::cout << " " << v[0] << " " << v[1] << " " << v[2] << " #";
        }
        std::cout << endl;
        
        ElementAccessor<3> ele_acc = ele.element_accessor();
        EXPECT_EQ(ele.dim(), ele_acc.dim());
        EXPECT_EQ(ele.centre()[0], ele_acc.centre()[0]);
        EXPECT_EQ(ele.centre()[1], ele_acc.centre()[1]);
        EXPECT_EQ(ele.centre()[2], ele_acc.centre()[2]);
    }
    
    delete mesh;
}








const string input = R"INPUT(
{   
   conc={ // formula on 3d 
       TYPE="FieldFormula",
       //value="log((x^2+y^2+z^2)^0.5)"
       value="x+y+z"
   },
   output_stream = {
    file = "test1.pvd", 
    format = {
        TYPE = "vtk", 
        variant = "ascii"
    }/*,
    output_mesh = {
        max_level = 3,
        refine_by_error = true,
        error_control_field = "conc"
    }*/
  }
  ,output = {fields = ["conc"]}
}
)INPUT";


TEST(OutputMesh, write_on_output_mesh) {
    
    typedef FieldAlgorithmBase<3, FieldValue<3>::Scalar > AlgScalarField;
    typedef Field<3,FieldValue<3>::Scalar> ScalarField;
  
    // setup FilePath directories
    FilePath::set_io_dirs(".",FilePath::get_absolute_working_dir(),"",".");

    // read mesh - simplset cube from test1
    FilePath mesh_file( "../mesh/simplest_cube.msh", FilePath::input_file);
    Mesh* mesh = new Mesh();
    ifstream in(string(mesh_file).c_str());
    mesh->read_gmsh_from_stream(in);
    
    
    // create scalar field out of FieldAlgorithmBase field
    ScalarField scalar_field;
    scalar_field.set_mesh(*mesh);
    
    // create field set of output fields
    EquationOutput output_fields;
    output_fields += scalar_field.name("conc").units(UnitSI::dimensionless()).flags_add(FieldFlag::allow_output);
    
    // create input record
    Input::Type::Record rec_type = Input::Type::Record("ErrorFieldTest","")
        .declare_key("conc", AlgScalarField::get_input_type_instance(), Input::Type::Default::obligatory(), "" )
        .declare_key("output_stream", OutputTime::get_input_type(), Input::Type::Default::obligatory(), "")
        .declare_key("output", output_fields.make_output_type("test_eq"), Input::Type::Default::obligatory(), "")
        .close();


    // read input string
    Input::ReaderToStorage reader( input, rec_type, Input::FileFormat::format_JSON );
    Input::Record in_rec=reader.get_root_interface<Input::Record>();

    // create FieldAlgorithmBase field
    auto alg_field = AlgScalarField::function_factory(in_rec.val<Input::AbstractRecord>("conc"), 1);
    
    // create field from FieldAlgorithmBase
    scalar_field.set_field(mesh->region_db().get_region_set("ALL"), alg_field, 0);
    scalar_field.output_type(OutputTime::NODE_DATA);
    
    // set time to all fields
    output_fields.set_time(0.0, LimitSide::right);
    
    // create output
    std::shared_ptr<OutputTime> output = std::make_shared<OutputVTK>();
    output->init_from_input("dummy_equation", *mesh, in_rec.val<Input::Record>("output_stream"));
    output_fields.initialize(output, in_rec.val<Input::Record>("output"), TimeGovernor());
    
    EXPECT_EQ(1,output_fields.size());
    EXPECT_TRUE(output_fields.is_field_output_time(*output_fields.field("conc"), 0.0));
    
    // register output fields, compute and write data
    output_fields.output(0.0);
    output->write_time_frame();

    delete mesh;
}









const string input_om = R"INPUT(
{   
    max_level = 2,
    refine_by_error = true,
    error_control_field = "conc"
}
)INPUT";

class TestOutputMesh : public testing::Test, public OutputMesh {
public:
    TestOutputMesh()
    : OutputMesh(mesh_, Input::ReaderToStorage( input_om, OutputMeshBase::get_input_type(), Input::FileFormat::format_JSON )
                            .get_root_interface<Input::Record>())
    {
    }
    
    ~TestOutputMesh()
    {
    }

    Mesh mesh_;
};


TEST_F(TestOutputMesh, read_input) {
    EXPECT_EQ(this->max_level_, 2);
    EXPECT_EQ(this->orig_mesh_, &(this->mesh_));
  
    Field<3,FieldValue<3>::Scalar> scalar_field;
    
    // create field set of output fields
    EquationOutput output_fields;
    output_fields += scalar_field.name("conc");
    this->select_error_control_field(output_fields);
    
    // no field 'conc' is to be found
    output_fields.field("conc")->name("conc_error");
    EXPECT_THROW( this->select_error_control_field(output_fields); ,
                  FieldSet::ExcUnknownField);
    
    // 'conc' field is now vector field
    Field<3,FieldValue<3>::VectorFixed> vector_field;
    output_fields += vector_field.name("conc");
    EXPECT_THROW( this->select_error_control_field(output_fields); ,
                  OutputMeshBase::ExcFieldNotScalar);
}
