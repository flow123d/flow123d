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
#include "fields/field.hh"
#include "fields/field_set.hh"
#include <io/equation_output.hh>
#include "fields/field_common.hh"
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
        xprintf(Msg,"%d %dD n_%d |",ele.idx(), ele.dim(), ele.n_nodes());
        for(unsigned int i=0; i < ele.n_nodes(); i++)
        {
            xprintf(Msg," %d",ele.node_index(i));
        }
        xprintf(Msg," |");
        for(auto& v : ele.vertex_list())
        {
            xprintf(Msg," %f %f %f #",v[0], v[1], v[2]);
        }
        xprintf(Msg,"\n");
        
        ElementAccessor<3> ele_acc = ele.element_accessor();
        EXPECT_EQ(ele.dim(), ele_acc.dim());
        EXPECT_EQ(ele.centre()[0], ele_acc.centre()[0]);
        EXPECT_EQ(ele.centre()[1], ele_acc.centre()[1]);
        EXPECT_EQ(ele.centre()[2], ele_acc.centre()[2]);
    }
    
    auto output_mesh_discont = std::make_shared<OutputMeshDiscontinuous>(*mesh);
    output_mesh_discont->create_mesh(output_mesh);
    
    xprintf(Msg,"DISCONTINUOUS\n");
    for(const auto &ele : *output_mesh_discont)
    {
        xprintf(Msg,"%d %dD n_%d |",ele.idx(), ele.dim(), ele.n_nodes());
        for(unsigned int i=0; i < ele.n_nodes(); i++)
        {
            xprintf(Msg," %d",ele.node_index(i));
        }
        xprintf(Msg," |");
        for(auto& v : ele.vertex_list())
        {
            xprintf(Msg," %f %f %f #",v[0], v[1], v[2]);
        }
        xprintf(Msg,"\n");
        
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
       value="log((x^2+y^2+z^2)^0.5)"
       //value="x+y+z"
   },
   output_stream = {
    file = "test1.pvd", 
    format = {
        TYPE = "vtk", 
        variant = "ascii"
    },
    output_mesh = {
        max_level = 3,
        refine_by_error = true,
        error_control_field = "conc"
    }
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
    //scalar_field.output_type(OutputTime::NODE_DATA);
    scalar_field.output_type(OutputTime::CORNER_DATA);
    
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

class TestOutputMesh : public testing::Test, public OutputMeshDiscontinuous {
public:
    TestOutputMesh()
    : OutputMeshDiscontinuous(mesh_, Input::ReaderToStorage( input_om, OutputMeshBase::get_input_type(), Input::FileFormat::format_JSON )
                            .get_root_interface<Input::Record>())
    {
    }
    
    template<int dim> void refine_single_element(AuxElement &el)
    {
        const int spacedim = 3;
        el.level = 0;
        std::vector<AuxElement> refs;
        
        std::vector<double> disc_coords;
        std::vector<unsigned int> disc_connectivity;
        
        this->refine_by_error_ = false;
        this->refine_aux_element<dim>(el, refs, ElementAccessor<3>(), nullptr);
        
        disc_coords.resize(spacedim * (dim+1) * refs.size());
        disc_connectivity.resize((dim+1) * refs.size());
        
        //gather coords and connectivity (disccontinous)
        for(unsigned int i=0; i < refs.size(); i++)
            for(unsigned int j=0; j < dim+1; j++)
            {
                unsigned int con = i*(dim+1) + j;
                disc_connectivity[con] = con;
                
                for(unsigned int k=0; k < spacedim; k++)
                    disc_coords[spacedim*con + k] = refs[i].nodes[j][k];
            }
            
        // correct data for the given aux element
        static const std::vector<double> res_coords[] = {
            {},
            { 0, 0, 0, 0.75, 0, 0, 0.75, 0, 0, 1.5, 0, 0, 1.5, 0, 0, 2.25, 0, 
            0, 2.25, 0, 0, 3, 0, 0 },
            {
            0, 0, 0, 0.75, 0, 0, 0, 0.75, 0, 0.75, 0, 0, 1.5, 0, 0, 0.75, 0.75, 
            0, 0, 0.75, 0, 0.75, 0.75, 0, 0, 1.5, 0, 0.75, 0, 0, 0.75, 0.75, 0, 
            0, 0.75, 0, 1.5, 0, 0, 2.25, 0, 0, 1.5, 0.75, 0, 2.25, 0, 0, 3, 
            0, 0, 2.25, 0.75, 0, 1.5, 0.75, 0, 2.25, 0.75, 0, 1.5, 1.5, 0, 2.25, 0, 
            0, 2.25, 0.75, 0, 1.5, 0.75, 0, 0, 1.5, 0, 0.75, 1.5, 0, 0, 2.25, 0, 
            0.75, 1.5, 0, 1.5, 1.5, 0, 0.75, 2.25, 0, 0, 2.25, 0, 0.75, 2.25, 0, 0, 
            3, 0, 0.75, 1.5, 0, 0.75, 2.25, 0, 0, 2.25, 0, 1.5, 0, 0, 1.5, 0.75, 
            0, 0.75, 0.75, 0, 1.5, 0.75, 0, 1.5, 1.5, 0, 0.75, 1.5, 0, 0.75, 0.75, 0, 
            0.75, 1.5, 0, 0, 1.5, 0, 1.5, 0.75, 0, 0.75, 1.5, 0, 0.75, 0.75, 0 },
            { 
            0, 0, 0, 0.75, 0, 0, 0, 0.75, 0, 0, 0, 0.75, 0.75, 0, 0, 1.5, 0, 
            0, 0.75, 0.75, 0, 0.75, 0, 0.75, 0, 0.75, 0, 0.75, 0.75, 0, 0, 1.5, 0, 
            0, 0.75, 0.75, 0, 0, 0.75, 0.75, 0, 0.75, 0, 0.75, 0.75, 0, 0, 1.5, 0.75, 
            0, 0, 0, 0, 0.75, 0.75, 0, 0.75, 0, 0.75, 0.75, 0.75, 0, 0, 0.75, 0.75, 
            0, 0.75, 0, 0.75, 0, 0.75, 0.75, 0.75, 0, 0, 0, 0, 0.75, 0, 0.75, 0, 
            0, 0.75, 0.75, 0.75, 0, 0, 0, 0.75, 0, 0.75, 0.75, 0, 0, 0.75, 0.75, 1.5, 
            0, 0, 2.25, 0, 0, 1.5, 0.75, 0, 1.5, 0, 0.75, 2.25, 0, 0, 3, 0, 
            0, 2.25, 0.75, 0, 2.25, 0, 0.75, 1.5, 0.75, 0, 2.25, 0.75, 0, 1.5, 1.5, 0, 
            1.5, 0.75, 0.75, 1.5, 0, 0.75, 2.25, 0, 0.75, 1.5, 0.75, 0.75, 1.5, 0, 1.5, 2.25, 
            0, 0, 1.5, 0, 0.75, 2.25, 0, 0.75, 1.5, 0.75, 0.75, 2.25, 0, 0, 2.25, 0.75, 
            0, 2.25, 0, 0.75, 1.5, 0.75, 0.75, 2.25, 0, 0, 1.5, 0, 0.75, 1.5, 0.75, 0, 
            1.5, 0.75, 0.75, 2.25, 0, 0, 1.5, 0.75, 0, 2.25, 0.75, 0, 1.5, 0.75, 0.75, 0, 
            1.5, 0, 0.75, 1.5, 0, 0, 2.25, 0, 0, 1.5, 0.75, 0.75, 1.5, 0, 1.5, 1.5, 
            0, 0.75, 2.25, 0, 0.75, 1.5, 0.75, 0, 2.25, 0, 0.75, 2.25, 0, 0, 3, 0, 
            0, 2.25, 0.75, 0, 1.5, 0.75, 0.75, 1.5, 0.75, 0, 2.25, 0.75, 0, 1.5, 1.5, 0.75, 
            1.5, 0, 0, 1.5, 0.75, 0.75, 1.5, 0.75, 0, 2.25, 0.75, 0.75, 1.5, 0, 0.75, 2.25, 
            0, 0.75, 1.5, 0.75, 0, 2.25, 0.75, 0.75, 1.5, 0, 0, 1.5, 0.75, 0, 2.25, 0, 
            0, 2.25, 0.75, 0.75, 1.5, 0, 0, 2.25, 0, 0.75, 2.25, 0, 0, 2.25, 0.75, 0, 
            0, 1.5, 0.75, 0, 1.5, 0, 0.75, 1.5, 0, 0, 2.25, 0.75, 0, 1.5, 1.5, 0, 
            1.5, 0.75, 0.75, 1.5, 0.75, 0, 2.25, 0, 0.75, 1.5, 0.75, 0.75, 1.5, 0, 1.5, 1.5, 
            0, 0.75, 2.25, 0, 0, 2.25, 0.75, 0, 2.25, 0, 0.75, 2.25, 0, 0, 3, 0.75, 
            0, 1.5, 0, 0, 2.25, 0.75, 0, 2.25, 0, 0.75, 2.25, 0.75, 0, 1.5, 0.75, 0.75, 
            1.5, 0.75, 0, 2.25, 0, 0.75, 2.25, 0.75, 0, 1.5, 0, 0, 2.25, 0, 0.75, 1.5, 
            0, 0.75, 2.25, 0.75, 0, 1.5, 0, 0.75, 1.5, 0.75, 0.75, 1.5, 0, 0.75, 2.25, 1.5, 
            0, 0, 0.75, 0, 0.75, 1.5, 0, 0.75, 0.75, 0.75, 0.75, 0.75, 0, 0.75, 0, 0, 
            1.5, 0.75, 0, 1.5, 0, 0.75, 1.5, 1.5, 0, 0.75, 0.75, 0, 1.5, 1.5, 0, 1.5, 
            0.75, 0.75, 1.5, 0.75, 0.75, 0.75, 0, 0.75, 1.5, 0.75, 0.75, 1.5, 0, 1.5, 1.5, 0.75, 
            0, 0.75, 0.75, 0.75, 0.75, 0, 0.75, 1.5, 0.75, 0.75, 1.5, 0.75, 0, 0.75, 0.75, 0, 
            1.5, 0, 0.75, 1.5, 0.75, 0.75, 1.5, 0.75, 0, 0.75, 0.75, 0.75, 0.75, 1.5, 0, 0.75, 
            0.75, 0.75, 1.5, 0.75, 0, 0.75, 1.5, 0, 0.75, 0.75, 0, 1.5, 0.75, 0.75, 1.5, 1.5, 
            0, 0, 1.5, 0.75, 0, 1.5, 0, 0.75, 0.75, 0.75, 0.75, 1.5, 0.75, 0, 1.5, 1.5, 
            0, 1.5, 0.75, 0.75, 0.75, 1.5, 0.75, 1.5, 0, 0.75, 1.5, 0.75, 0.75, 1.5, 0, 1.5, 
            0.75, 0.75, 1.5, 0.75, 0.75, 0.75, 0.75, 1.5, 0.75, 0.75, 0.75, 1.5, 0, 1.5, 1.5, 1.5, 
            0.75, 0, 0.75, 0.75, 0.75, 0.75, 1.5, 0.75, 0.75, 0.75, 1.5, 1.5, 0.75, 0, 1.5, 0.75, 
            0.75, 0.75, 1.5, 0.75, 0.75, 0.75, 1.5, 1.5, 0.75, 0, 0.75, 0.75, 0.75, 1.5, 0, 0.75, 
            0.75, 0.75, 1.5, 1.5, 0.75, 0, 1.5, 0, 0.75, 1.5, 0.75, 0.75, 0.75, 0.75, 1.5, 1.5, 
            0, 0, 0.75, 0, 0.75, 0.75, 0.75, 0, 0.75, 0.75, 0.75, 0.75, 0, 0.75, 0, 0, 
            1.5, 0, 0.75, 0.75, 0, 0.75, 1.5, 0.75, 0.75, 0, 0, 0.75, 0.75, 0, 1.5, 0, 
            0, 1.5, 0.75, 0.75, 0.75, 0.75, 0, 0.75, 1.5, 0, 1.5, 0.75, 0, 1.5, 1.5, 0.75, 
            0, 0.75, 0.75, 0.75, 0.75, 0, 0.75, 1.5, 0, 1.5, 0.75, 0.75, 0, 0.75, 0, 0.75, 
            0.75, 0, 0.75, 1.5, 0, 1.5, 0.75, 0.75, 0, 0.75, 0.75, 0.75, 0.75, 0.75, 0.75, 0, 
            0, 1.5, 0.75, 0.75, 0, 0.75, 0.75, 0.75, 0, 0, 0.75, 0.75, 0, 1.5, 0.75, 1.5, 
            0, 0, 0.75, 0.75, 0, 1.5, 0.75, 0, 0.75, 0.75, 0.75, 0.75, 0.75, 0, 0, 1.5, 
            0, 0.75, 1.5, 0, 0, 1.5, 0.75, 1.5, 0.75, 0, 0.75, 1.5, 0, 1.5, 1.5, 0, 
            0.75, 1.5, 0.75, 0.75, 0.75, 0.75, 0, 1.5, 0.75, 0.75, 1.5, 0.75, 0, 1.5, 1.5, 0.75, 
            0.75, 0, 0.75, 0.75, 0.75, 0, 1.5, 0.75, 0.75, 1.5, 0.75, 0.75, 0.75, 0, 0.75, 1.5, 
            0, 0, 1.5, 0.75, 0.75, 1.5, 0.75, 0.75, 0.75, 0, 0.75, 0.75, 0.75, 1.5, 0.75, 0, 
            0.75, 1.5, 0.75, 0.75, 0.75, 0, 1.5, 0.75, 0, 0.75, 1.5, 0, 0.75, 1.5, 0.75
            }
        };
        
        EXPECT_EQ(disc_coords.size(), res_coords[dim].size());
        
        for(unsigned int i=0; i < disc_coords.size(); i++)
        {
            EXPECT_DOUBLE_EQ(disc_coords[i], res_coords[dim][i]);
//             cout << disc_coords[i] << ", ";
//             if(i>0 && i%16 == 0) cout << endl;
        }
//         cout << "\n\n" << endl;
//         for(unsigned int i=0; i< disc_connectivity.size(); i++)
//         {
//             cout << i << "  "; 
//             for(unsigned int k=0; k<spacedim; k++){
//                 cout << disc_coords[i*spacedim+k] << " ";
//             }
//             cout << endl;
//         } 
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
}




TEST_F(TestOutputMesh, refine) {
    
    const double a = 3.0;
    
    AuxElement el1;
    el1.nodes = {{0, 0, 0}, {a, 0, 0}};
    this->refine_single_element<1>(el1);
    
    AuxElement el2;
    el2.nodes = {{0, 0, 0}, {a, 0, 0}, {0, a, 0}};
    this->refine_single_element<2>(el2);
   
    AuxElement el3;
    el3.nodes = {{0, 0, 0}, {a, 0, 0}, {0, a, 0}, {0, 0, a}};
    this->refine_single_element<3>(el3);
    
}



const string input_refine = R"INPUT(
{   
   conc={ // formula on 3d 
       TYPE="FieldFormula",
       value="if((x^2+y^2+z^2)^0.5 > 0.01, log((x^2+y^2+z^2)^0.5), log(0.01))"
   },
   output_stream = {
    file = "test_refine.pvd",
    format = {
        TYPE = "vtk", 
        variant = "ascii"
    },
    output_mesh = {
        max_level = 3
        refine_by_error = true
        error_control_field = "conc"
    }
  }
  ,output = {fields = ["conc"]}
}
)INPUT";

// simplest mesh
string small_mesh = R"CODE(
$MeshFormat
2.2 0 8
$EndMeshFormat
$Nodes
4
1 -3 -2 0
2 2 -3 0
3 3 2 0
4 -2 3 0
$EndNodes
$Elements
2
1 2 2 39 40 1 2 3
2 2 2 39 40 1 3 4
$EndElements
)CODE";

TEST(OutputMesh, write_on_refined_mesh) {
    
    typedef FieldAlgorithmBase<3, FieldValue<3>::Scalar > AlgScalarField;
    typedef Field<3,FieldValue<3>::Scalar> ScalarField;
  
    // setup FilePath directories
    FilePath::set_io_dirs(".",FilePath::get_absolute_working_dir(),"",".");

    // read mesh - simplset cube from test1
    Mesh* mesh = new Mesh();
    stringstream in(small_mesh.c_str());
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
    Input::ReaderToStorage reader( input_refine, rec_type, Input::FileFormat::format_JSON );
    Input::Record in_rec=reader.get_root_interface<Input::Record>();

    // create FieldAlgorithmBase field
    auto alg_field = AlgScalarField::function_factory(in_rec.val<Input::AbstractRecord>("conc"), 1);
    
    // create field from FieldAlgorithmBase
    scalar_field.set_field(mesh->region_db().get_region_set("ALL"), alg_field, 0);
//     scalar_field.output_type(OutputTime::NODE_DATA);
    scalar_field.output_type(OutputTime::CORNER_DATA);
    
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


// class TestOutputMesh : public testing::Test, public OutputMesh {
// public:
//     TestOutputMesh()
//     : OutputMesh(nullptr,2)
//     {
//     }
//     
//     template<int dim> void refine_single_element(AuxElement &el)
//     {
//         const int spacedim = 3;
//         el.level = 0;
//         unsigned int last = dim;
//         std::vector<AuxElement> refs;
//         
//         std::vector<double> coords;
//         std::vector<unsigned int> connectivity;
//         
//         this->refine_aux_element<dim>(el, refs, last);
//         
//         coords.resize(spacedim * last, std::numeric_limits<double>::quiet_NaN());
//         connectivity.resize(refs.size()*(dim+1));
//         
//         //gather coords and connectivity (in a continous way inside element)
//         for(unsigned int i=0; i < refs.size(); i++)
//             for(unsigned int j=0; j < dim+1; j++)
//             {
//                 unsigned int con = refs[i].connectivity[j];
//                 connectivity[(dim+1)*i + j] = con;
//                 
//                 if(std::isnan(coords[spacedim*con])) //if we haven't add this coords
//                 for(unsigned int k=0; k < spacedim; k++)
//                     coords[spacedim*con + k] = refs[i].nodes[j][k];
//             }
//             
//         // correct data for the given aux element
//         static const std::vector<unsigned int> res_conn[] = {
//             {},
//             { 0, 3, 3, 2, 2, 4, 4, 1 },
//             {
//             0, 6, 7, 6, 3, 8, 7, 8, 4, 6, 8, 7, 3, 
//             9, 10, 9, 1, 11, 10, 11, 5, 9, 11, 10, 4, 
//             12, 13, 12, 5, 14, 13, 14, 2, 12, 14, 13, 3, 
//             15, 16, 15, 5, 17, 16, 17, 4, 15, 17, 16 },
//             {
//             0, 10, 11, 3, 10, 4, 12, 14, 11, 12, 5, 15, 13, 14, 15, 3, 10, 
//             13, 14, 15, 10, 12, 14, 15, 10, 13, 11, 15, 10, 11, 12, 15, 4, 
//             16, 17, 8, 16, 1, 18, 20, 17, 18, 6, 21, 19, 20, 21, 8, 16, 
//             19, 20, 21, 16, 18, 20, 21, 16, 19, 17, 21, 16, 17, 18, 21, 5, 
//             22, 23, 9, 22, 6, 24, 26, 23, 24, 2, 27, 25, 26, 27, 9, 22, 
//             25, 26, 27, 22, 24, 26, 27, 22, 25, 23, 27, 22, 23, 24, 27, 7, 
//             28, 29, 3, 28, 8, 30, 32, 29, 30, 9, 33, 31, 32, 33, 3, 28, 
//             31, 32, 33, 28, 30, 32, 33, 28, 31, 29, 33, 28, 29, 30, 33, 4, 
//             34, 35, 9, 34, 7, 36, 38, 35, 36, 8, 39, 37, 38, 39, 9, 34, 
//             37, 38, 39, 34, 36, 38, 39, 34, 37, 35, 39, 34, 35, 36, 39, 4, 
//             40, 41, 9, 40, 6, 42, 44, 41, 42, 8, 45, 43, 44, 45, 9, 40, 
//             43, 44, 45, 40, 42, 44, 45, 40, 43, 41, 45, 40, 41, 42, 45, 4, 
//             46, 47, 9, 46, 7, 48, 50, 47, 48, 5, 51, 49, 50, 51, 9, 46, 
//             49, 50, 51, 46, 48, 50, 51, 46, 49, 47, 51, 46, 47, 48, 51, 4, 
//             52, 53, 9, 52, 5, 54, 56, 53, 54, 6, 57, 55, 56, 57, 9, 52, 
//             55, 56, 57, 52, 54, 56, 57, 52, 55, 53, 57, 52, 53, 54, 57 }
//         };
//         static const std::vector<double> res_coords[] = {
//             {},
//             { 0, 0, 0, 3, 0, 0, 1.5, 0, 0, 0.75, 0, 0 },
//             {
//             0, 0, 0, 3, 0, 0, 0, 3, 0, 1.5, 
//             0, 0, 0, 1.5, 0, 1.5, 1.5, 0, 0.75, 
//             0, 0, 0, 0.75, 0, 0.75, 0.75, 0, 2.25, 
//             0, 0, 1.5, 0.75, 0, 2.25, 0.75, 0, 0.75, 
//             1.5, 0, 0, 2.25, 0, 0.75, 2.25, 0, 1.5, 
//             0.75, 0, 0.75, 0.75, 0 },
//             {
//             0, 0, 0, 3, 0, 0, 0, 3, 0, 0, 0, 3, 1.5, 0, 0, 0, 1.5, 
//             0, 1.5, 1.5, 0, 0, 0, 1.5, 1.5, 0, 1.5, 0, 1.5, 1.5, 0.75, 0, 0, 
//             0, 0.75, 0, 0.75, 0.75, 0, 0, 0, 1.5, 0.75, 0, 1.5, 0, 0.75, 1.5, 2.25, 
//             0, 0, 1.5, 0.75, 0, 2.25, 0.75, 0, 1.5, 0, 0.75, 2.25, 0, 0.75, 1.5, 0.75, 
//             0.75, 0.75, 1.5, 0, 0, 2.25, 0, 0.75, 2.25, 0, 0, 1.5, 0.75, 0.75, 1.5, 0.75, 
//             0, 2.25, 0.75, 0.75, 0, 1.5, 0, 0.75, 1.5, 0.75, 0.75, 1.5, 0, 0, 2.25, 0.75, 
//             0, 2.25, 0, 0.75, 2.25, 0.75, 0, 0.75, 1.5, 0, 0.75, 0.75, 0, 1.5, 0.75, 0.75, 
//             0.75, 0, 0.75, 1.5, 0.75, 0.75, 1.5, 1.5, 0.75, 0, 1.5, 0, 0.75, 1.5, 0.75, 0.75, 
//             0.75, 0.75, 0.75, 0.75, 1.5, 0.75, 0.75, 0.75, 1.5, 0.75, 0, 0.75, 0.75, 0.75, 0, 0, 
//             0.75, 0.75, 0.75, 0.75, 0.75, 0, 0.75, 1.5, 0, 1.5, 0.75, 0.75, 0.75, 0, 1.5, 0.75, 
//             0, 0.75, 1.5, 0, 0.75, 0.75, 0.75, 0, 1.5, 0.75 }
//         };
//         
//         EXPECT_EQ(coords.size(), res_coords[dim].size());
//         EXPECT_EQ(connectivity.size(), res_conn[dim].size());
//         for(unsigned int i=0; i < connectivity.size(); i++)
//         {
//             EXPECT_EQ(connectivity[i], res_conn[dim][i]);
// //             cout << connectivity[i] << ", ";
// //             if(i>0 && i%16 == 0) cout << endl;
//         }
// //         cout << "\n\n" <<  endl;
//         for(unsigned int i=0; i < coords.size(); i++)
//         {
//             EXPECT_DOUBLE_EQ(coords[i], res_coords[dim][i]);
// //             cout << coords[i] << ", ";
// //             if(i>0 && i%16 == 0) cout << endl;
//         }
// //         cout << endl;
//         
//         // just text output
// //         for(unsigned int i=0; i < refs.size(); i++)
// //         {   
// //             for(unsigned int j=0; j < dim+1; j++)
// //                 cout << refs[i].connectivity[j] << " ";
// //             cout << endl;    
// //             
// //             for(unsigned int j=0; j < dim+1; j++)
// //             {
// //                 refs[i].nodes[j].print();
// //                 cout << "-------------------------------\n";
// //             }
// //             cout << "####################################\n";
// //         }
//     }
// 
//     ~TestOutputMesh()
//     {
//     }
// };
// 
// 
// TEST_F(TestOutputMesh, refine) {
//     
//     const double a = 3.0;
//     
//     AuxElement el1;
//     el1.nodes = {{0, 0, 0}, {a, 0, 0}};
//     el1.connectivity = {0, 1};
//     this->refine_single_element<1>(el1);
//     
//     AuxElement el2;
//     el2.nodes = {{0, 0, 0}, {a, 0, 0}, {0, a, 0}};
//     el2.connectivity = {0, 1, 2};
//     this->refine_single_element<2>(el2);
//    
//     AuxElement el3;
//     el3.nodes = {{0, 0, 0}, {a, 0, 0}, {0, a, 0}, {0, 0, a}};
//     el3.connectivity = {0, 1, 2, 3};
//     this->refine_single_element<3>(el3);
//     
// }