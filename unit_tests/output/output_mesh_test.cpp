/*
 * output_vtk_test.cpp
 *
 *  Created on: June 9, 2016
 *      Author: pe
 */

#define TEST_USE_PETSC
#define FEAL_OVERRIDE_ASSERTS
#include <flow_gtest_mpi.hh>
#include <mesh_constructor.hh>
#include "io/output_time.hh"
#include "io/output_vtk.hh"
#include "mesh/mesh.h"
#include "io/msh_gmshreader.h"
#include "input/reader_to_storage.hh"
#include "system/sys_profiler.hh"

#include "fields/field_constant.hh"
#include "fields/field.hh"
#include "fields/field_set.hh"
#include <fields/equation_output.hh>
#include "fields/field_common.hh"
#include "input/input_type.hh"

#include "input/accessors.hh"

#include "io/output_mesh.hh"
#include "io/output_element.hh"


FLOW123D_FORCE_LINK_IN_PARENT(field_formula)



TEST(OutputMesh, create_identical)
{
    // read mesh - simplest cube from test1
    FilePath mesh_file( string(UNIT_TESTS_SRC_DIR) + "/mesh/simplest_cube.msh", FilePath::input_file);
    Mesh *mesh = mesh_full_constructor("{mesh_file=\"" + (string)mesh_file + "\"}");
    
    auto output_mesh = std::make_shared<OutputMesh>(*mesh);
    output_mesh->create_identical_mesh();
    
    std::cout << "nodes: ";
    output_mesh->nodes_->print_ascii_all(std::cout);
    std::cout << endl;
    std::cout << "connectivity: ";
    output_mesh->connectivity_->print_ascii_all(std::cout);
    std::cout << endl;
    std::cout << "offsets: ";
    output_mesh->offsets_->print_ascii_all(std::cout);
    std::cout << endl;
    
    EXPECT_EQ(output_mesh->nodes_->n_values(), mesh->n_nodes());
    EXPECT_EQ(output_mesh->offsets_->n_values(), mesh->n_elements());
    
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
       value="1/((x^2+y^2+z^2)^0.5)"
       //value="log((x^2+y^2+z^2)^0.5)"
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
    Mesh *mesh = mesh_full_constructor("{mesh_file=\"../mesh/simplest_cube.msh\"}");
    
    
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
    FieldAlgoBaseInitData init_data("conc", 1, UnitSI::dimensionless());
    auto alg_field = AlgScalarField::function_factory(in_rec.val<Input::AbstractRecord>("conc"), init_data);
    
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
    : OutputMeshDiscontinuous(*mesh_, Input::ReaderToStorage( input_om, 
                                                              const_cast<Input::Type::Record &>(
                                                                    OutputMeshBase::get_input_type()),
                                                              Input::FileFormat::format_JSON )
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
            1.5, 1.5, 0, 1.5, 0.75, 0, 2.25, 0.75, 0, 1.5, 0.75, 0.75, 1.5, 0.75, 0, 1.5, 0, 
            0, 2.25, 0, 0, 1.5, 0, 0.75, 2.25, 0.75, 0, 2.25, 0, 0, 3, 0, 0, 
            2.25, 0, 0.75, 1.5, 0.75, 0.75, 1.5, 0, 0.75, 2.25, 0, 0.75, 1.5, 0, 1.5, 1.5, 
            0.75, 0, 2.25, 0.75, 0, 1.5, 0.75, 0.75, 1.5, 0, 0.75, 1.5, 0.75, 0, 2.25, 0.75, 
            0, 2.25, 0, 0, 1.5, 0, 0.75, 2.25, 0.75, 0, 1.5, 0.75, 0.75, 1.5, 0, 0.75, 
            2.25, 0, 0.75, 2.25, 0.75, 0, 2.25, 0, 0, 1.5, 0, 0.75, 2.25, 0, 0.75, 0, 
            3, 0, 0, 2.25, 0, 0.75, 2.25, 0, 0, 2.25, 0.75, 0, 2.25, 0, 0, 1.5, 
            0, 0.75, 1.5, 0, 0, 1.5, 0.75, 0.75, 2.25, 0, 0.75, 1.5, 0, 1.5, 1.5, 0, 
            0.75, 1.5, 0.75, 0, 2.25, 0.75, 0, 1.5, 0.75, 0.75, 1.5, 0.75, 0, 1.5, 1.5, 0, 
            2.25, 0, 0.75, 2.25, 0, 0, 2.25, 0.75, 0, 1.5, 0.75, 0, 2.25, 0, 0.75, 2.25, 
            0, 0.75, 1.5, 0, 0, 1.5, 0.75, 0.75, 2.25, 0, 0, 2.25, 0.75, 0, 1.5, 0.75, 
            0.75, 1.5, 0.75, 0.75, 2.25, 0, 0.75, 1.5, 0, 0, 1.5, 0.75, 0.75, 1.5, 0.75, 0, 
            1.5, 0, 0, 0.75, 0, 0.75, 0.75, 0, 0, 0.75, 0.75, 0, 0.75, 0, 0, 0, 
            0, 0.75, 0, 0, 0, 0, 0.75, 0.75, 0.75, 0, 0.75, 0, 0, 1.5, 0, 0, 
            0.75, 0, 0.75, 0, 0.75, 0.75, 0, 0, 0.75, 0.75, 0, 0.75, 0, 0, 1.5, 0, 
            0.75, 0, 0.75, 0.75, 0, 0, 0.75, 0.75, 0, 0, 0.75, 0, 0.75, 0, 0.75, 0.75, 
            0, 0.75, 0, 0, 0, 0, 0.75, 0.75, 0.75, 0, 0, 0.75, 0.75, 0, 0, 0.75, 
            0.75, 0, 0.75, 0.75, 0.75, 0, 0.75, 0, 0, 0, 0, 0.75, 0.75, 0, 0.75, 0, 
            1.5, 1.5, 0, 0.75, 1.5, 0.75, 0.75, 1.5, 0, 0.75, 2.25, 0, 0.75, 1.5, 0, 0, 
            1.5, 0.75, 0, 1.5, 0, 0, 2.25, 0.75, 0.75, 1.5, 0.75, 0, 1.5, 1.5, 0, 1.5, 
            0.75, 0, 2.25, 0, 0.75, 2.25, 0, 0, 2.25, 0.75, 0, 2.25, 0, 0, 3, 0, 
            0.75, 1.5, 0.75, 0.75, 1.5, 0, 0.75, 2.25, 0, 0, 2.25, 0, 0.75, 1.5, 0.75, 0.75, 
            1.5, 0.75, 0, 1.5, 0, 0, 2.25, 0.75, 0.75, 1.5, 0, 0.75, 2.25, 0, 0, 2.25, 
            0.75, 0, 2.25, 0.75, 0.75, 1.5, 0.75, 0, 1.5, 0, 0, 2.25, 0.75, 0, 2.25, 1.5, 
            0, 0, 1.5, 0, 0.75, 1.5, 0.75, 0, 0.75, 0.75, 0.75, 1.5, 0, 0.75, 1.5, 0, 
            1.5, 1.5, 0.75, 0.75, 0.75, 0.75, 1.5, 1.5, 0.75, 0, 1.5, 0.75, 0.75, 1.5, 1.5, 0, 
            0.75, 1.5, 0.75, 0.75, 0.75, 0.75, 0.75, 0.75, 1.5, 0.75, 1.5, 0.75, 0, 1.5, 1.5, 0.75, 
            0.75, 0.75, 1.5, 0.75, 0.75, 1.5, 0.75, 0, 0.75, 1.5, 0.75, 1.5, 0.75, 0.75, 0.75, 0.75, 
            0.75, 0.75, 0.75, 1.5, 0.75, 1.5, 0.75, 1.5, 0.75, 0.75, 0.75, 0.75, 0.75, 1.5, 0.75, 0, 
            1.5, 0, 0.75, 0.75, 0.75, 0.75, 1.5, 0.75, 0.75, 0.75, 0.75, 1.5, 1.5, 0, 0.75, 1.5, 
            0, 0, 0.75, 0.75, 0, 1.5, 0.75, 0, 0.75, 0.75, 0.75, 0.75, 0.75, 0, 0, 1.5, 
            0, 0.75, 1.5, 0, 0, 1.5, 0.75, 1.5, 0.75, 0, 0.75, 1.5, 0, 1.5, 1.5, 0, 
            0.75, 1.5, 0.75, 0.75, 0.75, 0.75, 0, 1.5, 0.75, 0.75, 1.5, 0.75, 0, 1.5, 1.5, 0.75, 
            0.75, 0.75, 0.75, 1.5, 0, 1.5, 0.75, 0, 0.75, 1.5, 0.75, 0.75, 1.5, 0, 0.75, 0.75, 
            0.75, 0, 1.5, 0.75, 0.75, 1.5, 0.75, 0.75, 1.5, 0, 0.75, 0.75, 0.75, 1.5, 0.75, 0, 
            0.75, 0.75, 0, 0.75, 0.75, 0.75, 0.75, 1.5, 0, 0, 1.5, 0.75, 0.75, 0.75, 0, 1.5, 
            0, 1.5, 0.75, 0.75, 1.5, 1.5, 0, 0.75, 0.75, 0, 1.5, 0.75, 0.75, 1.5, 0, 1.5, 
            1.5, 0.75, 0.75, 0.75, 0, 0.75, 1.5, 1.5, 0, 0.75, 0.75, 0.75, 0.75, 1.5, 0, 0, 
            0.75, 0, 0.75, 0.75, 0, 1.5, 0, 0.75, 1.5, 0.75, 0, 0.75, 0, 0, 1.5, 0.75, 
            0, 1.5, 0.75, 0.75, 0.75, 1.5, 0, 0.75, 0.75, 0, 0.75, 0.75, 0.75, 0.75, 0.75, 0, 
            1.5, 0, 0.75, 1.5, 0.75, 0, 0.75, 0.75, 0.75, 0.75, 0.75, 0, 1.5, 1.5, 0, 0.75, 
            0.75, 0.75, 1.5, 0.75, 0, 1.5, 0.75, 0.75, 0.75, 0, 0.75, 1.5, 0.75, 0.75, 1.5, 0, 
            1.5, 0, 0, 1.5, 0.75, 0.75, 0.75, 0, 0, 0.75, 0.75, 0, 1.5, 0.75, 0, 1.5, 
            1.5, 0.75, 0.75, 0.75, 0, 0.75, 1.5, 0.75, 0.75, 0, 0.75, 0.75, 0.75, 1.5, 0, 0, 
            0.75, 0, 0.75, 0, 0.75, 0.75, 0, 0.75, 1.5, 0.75, 0, 0.75, 0, 0, 1.5, 0, 
            0.75, 0.75, 0.75, 0.75, 0.75, 0.75, 0.75, 0, 0.75, 0, 0.75, 0.75, 0.75, 0.75, 0, 0.75, 
            0.75, 0, 0.75, 1.5, 0.75, 0, 0.75, 0.75, 0.75, 0.75, 0, 0.75, 0.75, 0.75, 0.75, 0, 
            0, 1.5, 0.75, 0, 0.75, 0.75, 0.75, 0.75, 0.75, 0, 0.75, 1.5, 0, 1.5, 0.75
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

    Mesh * mesh_;
};



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

TEST(OutputMesh, write_on_refined_mesh) {
    
    typedef FieldAlgorithmBase<3, FieldValue<3>::Scalar > AlgScalarField;
    typedef Field<3,FieldValue<3>::Scalar> ScalarField;
    
    // setup FilePath directories
    FilePath::set_io_dirs(".",FilePath::get_absolute_working_dir(),"",".");

    // read mesh - simplset cube from test1
    Mesh *mesh = mesh_full_constructor("{mesh_file=\"../mesh/simplest_cube.msh\"}");
    
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
    FieldAlgoBaseInitData field_algo_init_data("conc",1, UnitSI::one());
    auto alg_field = AlgScalarField::function_factory(in_rec.val<Input::AbstractRecord>("conc"), field_algo_init_data);
    
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
    
    //test select_error_control_field exception
    // no field 'conc' is to be found
    output_fields.field("conc")->name("conc_error");
    EXPECT_THROW( output_fields.select_error_control_field(); ,
                  FieldSet::ExcUnknownField);
    output_fields.field("conc_error")->name("conc");
    
    // register output fields, compute and write data
    output_fields.output(0.0);
    output->write_time_frame();

    delete mesh;
}
