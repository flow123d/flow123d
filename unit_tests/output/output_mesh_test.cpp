/*
 * output_vtk_test.cpp
 *
 *  Created on: June 9, 2016
 *      Author: pe
 */

#define TEST_USE_PETSC
#define FEAL_OVERRIDE_ASSERTS
#include <flow_gtest_mpi.hh>
#include "arma_expect.hh"

#include <mesh_constructor.hh>
#include "mesh/mesh.h"
#include "input/input_type.hh"
#include "input/accessors.hh"

#include "io/output_mesh.hh"
#include "io/output_element.hh"


class OutputMeshTest : public OutputMesh {
public:
	OutputMeshTest(Mesh  &mesh)
	: OutputMesh(mesh) {}

    void print_and_check_sizes(Mesh *mesh)
    {
        std::cout << "nodes: ";
        this->nodes_->print_ascii_all(std::cout);
        std::cout << endl;
        std::cout << "connectivity: ";
        this->connectivity_->print_ascii_all(std::cout);
        std::cout << endl;
        std::cout << "offsets: ";
        this->offsets_->print_ascii_all(std::cout);
        std::cout << endl;

        EXPECT_EQ(this->nodes_->n_values(), mesh->n_nodes());
        EXPECT_EQ(this->offsets_->n_values(), mesh->n_elements());
    }

    ~OutputMeshTest() {}

};


TEST(OutputMesh, create_identical)
{
    // read mesh - simplest cube from test1
    FilePath mesh_file( string(UNIT_TESTS_SRC_DIR) + "/mesh/simplest_cube.msh", FilePath::input_file);
    Mesh *mesh = mesh_full_constructor("{mesh_file=\"" + (string)mesh_file + "\"}");
    
    auto output_mesh = std::make_shared<OutputMeshTest>(*mesh);
    output_mesh->create_sub_mesh();
    output_mesh->make_serial_master_mesh();
    output_mesh->print_and_check_sizes(mesh);
    
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
        EXPECT_EQ(ele_acc.dim(), ele.dim());
        EXPECT_ARMA_EQ(ele_acc.centre(), ele.centre());
    }
    
    auto output_mesh_discont = std::make_shared<OutputMeshDiscontinuous>(*mesh);
    output_mesh_discont->create_sub_mesh();
    output_mesh_discont->make_serial_master_mesh();
    
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
        EXPECT_EQ(ele_acc.dim(), ele.dim());
        EXPECT_ARMA_EQ(ele_acc.centre(), ele.centre());
    }
    
    delete mesh;
}



const string input_om = R"INPUT(
{   
    max_level = 2
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
        this->refine_aux_element<dim>(el, refs, ElementAccessor<3>());
        
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
