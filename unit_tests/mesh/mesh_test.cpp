/*
 * mesh_test.cpp
 *
 *  Created on: Jan 26, 2013
 *      Author: jb
 */

#define TEST_USE_PETSC
#define FEAL_OVERRIDE_ASSERTS
#include <flow_gtest_mpi.hh>
#include <mesh_constructor.hh>

#include "mesh/mesh.h"
#include "mesh/bc_mesh.hh"
#include "io/msh_gmshreader.h"
#include <iostream>
#include <vector>
#include "mesh/accessors.hh"
#include "mesh/partitioning.hh"
#include "input/reader_to_storage.hh"
#include "system/sys_profiler.hh"



using namespace std;

class MeshTest :  public testing::Test, public Mesh {
public:
    MeshTest()
    : Mesh()
    {
    }

    ~MeshTest()
    {
    }
};


TEST_F(MeshTest, intersect_nodes_lists) {
	node_elements_.resize(3);
	node_elements_[0]={ 0, 1, 2, 3, 4};
	node_elements_[1]={ 0, 2, 3, 4};
	node_elements_[2]={ 0, 1, 2, 4};

    vector<unsigned int> node_list={0,1,2};
    vector<unsigned int> result;
    intersect_element_lists(node_list, result);
    EXPECT_EQ( vector<unsigned int>( {0,2,4} ), result );

    node_list={0,1};
    intersect_element_lists(node_list, result);
    EXPECT_EQ( vector<unsigned int>( {0,2,3,4} ), result );

    node_list={0};
    intersect_element_lists(node_list, result);
    EXPECT_EQ( vector<unsigned int>( {0,1,2,3,4} ), result );

}



TEST(MeshTopology, make_neighbours_and_edges) {
	// has to introduce some flag for passing absolute path to 'test_units' in source tree
	FilePath::set_io_dirs(".",UNIT_TESTS_SRC_DIR,"",".");

    Profiler::instance();
    
    Mesh * mesh = mesh_full_constructor("{mesh_file=\"mesh/simplest_cube.msh\"}");

    EXPECT_EQ(9, mesh->n_elements());
    EXPECT_EQ(18, mesh->n_elements(true));

    // check boundary elements
    EXPECT_EQ(101 , mesh->element_accessor(9).region().id() );
    EXPECT_EQ(101 , mesh->element_accessor(10).region().id() );
    EXPECT_EQ(102 , mesh->element_accessor(11).region().id() );
    EXPECT_EQ(102 , mesh->element_accessor(12).region().id() );
    EXPECT_EQ( -3 , int( mesh->element_accessor(13).region().id() ) );
    EXPECT_EQ( -3 , int( mesh->element_accessor(26).region().id() ) );

    //check edges
    EXPECT_EQ(28,mesh->n_edges());

    //check neighbours
    EXPECT_EQ(6, mesh->n_vb_neighbours() );

    delete mesh;
}


const string mesh_input = R"YAML(
mesh_file: "mesh/simplest_cube.msh"
regions:
 - !From_Elements
   id: 3000
   name: new region
   element_list:
    - 6
    - 7
 - !From_Id
   id: 37
   name: 1D rename
 - !Union
   name: 3D
   region_ids:
    - 39
    - 40
)YAML";

TEST(Mesh, init_from_input) {
    FilePath::set_io_dirs(".",UNIT_TESTS_SRC_DIR,"",".");

    Mesh * mesh = mesh_full_constructor(mesh_input, Input::FileFormat::format_YAML);

    EXPECT_EQ( 37, mesh->element_accessor(0).region().id() );
    EXPECT_EQ( "1D rename", mesh->element_accessor(0).region().label() );

    EXPECT_EQ( 38, mesh->element_accessor(1).region().id() );
    EXPECT_EQ( 38, mesh->element_accessor(2).region().id() );
    EXPECT_EQ( 39, mesh->element_accessor(3).region().id() );
    EXPECT_EQ( 39, mesh->element_accessor(4).region().id() );
    EXPECT_EQ( 3000, mesh->element_accessor(5).region().id() );
    EXPECT_EQ( 3000, mesh->element_accessor(6).region().id() );
    EXPECT_EQ( 40, mesh->element_accessor(7).region().id() );
    EXPECT_EQ( 40, mesh->element_accessor(8).region().id() );

    RegionSet set = mesh->region_db().get_region_set("3D");
    EXPECT_EQ( 39, set[0].id() );
    EXPECT_EQ( 40, set[1].id() );

    delete mesh;
}


TEST(Mesh, decompose_problem) {
	FilePath::set_io_dirs(".",UNIT_TESTS_SRC_DIR,"",".");

	std::string mesh_in_string = "{mesh_file=\"mesh/decompose_problem.msh\"}";
	Mesh * mesh = mesh_constructor(mesh_in_string);
    auto reader = reader_constructor(mesh_in_string);
    reader->read_physical_names(mesh);
    reader->read_raw_mesh(mesh);

    int np;
    MPI_Comm_size(MPI_COMM_WORLD, &np);
    if (np > 1) {
        EXPECT_THROW_WHAT( { mesh->setup_topology(); }, Partitioning::ExcDecomposeMesh,
                "greater then number of elements 1. Can not make partitioning of the mesh");
    }

    delete mesh;
}


TEST(Mesh, check_compatible_mesh) {
	FilePath::set_io_dirs(".",UNIT_TESTS_SRC_DIR,"",".");

	vector<LongIdx> bulk_elms_id, boundary_elms_id;

    std::string mesh_string = "{mesh_file=\"mesh/simplest_cube.msh\"}";
    Mesh * target_mesh = mesh_constructor(mesh_string);
    auto target_reader = reader_constructor(mesh_string);
    target_reader->read_physical_names(target_mesh);
    target_reader->read_raw_mesh(target_mesh);
    target_mesh->setup_topology();

    {
        std::string mesh_in_string = "{mesh_file=\"mesh/pvd-test/pvd-test-000000.vtu\"}";
        Mesh * mesh = mesh_constructor(mesh_in_string);
        auto reader = reader_constructor(mesh_in_string);
        //reader->read_physical_names(mesh); // not implemented
        reader->read_raw_mesh(mesh);

        EXPECT_TRUE( mesh->check_compatible_mesh(*target_mesh, bulk_elms_id, boundary_elms_id) );

        delete mesh;
    }

    {
        std::string mesh_in_string = "{mesh_file=\"mesh/test_108_elem.msh\"}";
        Mesh * mesh = mesh_constructor(mesh_in_string);
        auto reader = reader_constructor(mesh_in_string);
        // reader->read_physical_names(mesh); // not implemented
        reader->read_raw_mesh(mesh);

        EXPECT_FALSE( mesh->check_compatible_mesh(*target_mesh, bulk_elms_id, boundary_elms_id) );

        delete mesh;
    }

    delete target_mesh;
}


TEST(BCMesh, element_ranges) {
	FilePath::set_io_dirs(".",UNIT_TESTS_SRC_DIR,"",".");

	std::string mesh_in_string = "{mesh_file=\"mesh/simplest_cube.msh\"}";
	Mesh * mesh = mesh_constructor(mesh_in_string);
    auto reader = reader_constructor(mesh_in_string);
    reader->read_physical_names(mesh);
    reader->read_raw_mesh(mesh);

    BCMesh *bc_mesh = mesh->get_bc_mesh();
    unsigned int expected_val = 0;

    for (auto elm : mesh->elements_range()) {
    	EXPECT_EQ(elm.idx(), expected_val);
    	EXPECT_EQ(elm.mesh_idx(), expected_val);
    	expected_val++;
    }
    for (auto elm : bc_mesh->elements_range()) {
    	EXPECT_EQ(elm.idx(), expected_val-mesh->n_elements());
    	EXPECT_EQ(elm.mesh_idx(), expected_val);
    	expected_val++;
    }

    //delete bc_mesh;
    delete mesh;
}
