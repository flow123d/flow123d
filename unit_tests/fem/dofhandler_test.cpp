#define TEST_USE_PETSC
#include <flow_gtest_mpi.hh>
#include <cmath>
#include "fem/fe_p.hh"
#include "mesh/mesh.h"
#include <mesh_constructor.hh>
#include "fem/dofhandler.hh"


// simple mesh consisting of 4 triangles and 1 line element
//
//  4---------------5
//  | \           / |
//  |   \   5   1   |
//  |     \   /     |
//  |  3    3    4  |
//  |     /   \     |
//  |   /   2   \   |
//  | /           \ |
//  1---------------2
string small_mesh = R"CODE(
$MeshFormat
2.2 0 8
$EndMeshFormat
$Nodes
5
1 0 0 0
2 2 0 0
3 1 1 0
4 0 2 0
5 2 2 0
$EndNodes
$Elements
5
1 1 2 1 1 3 5
2 2 2 2 2 1 2 3
3 2 2 2 2 1 3 4
4 2 2 2 2 3 2 5
5 2 2 2 2 3 5 4
$EndElements
)CODE";


// simple mesh consisting of 4 triangles and 1 line element
//
//  5---------------6
//  | \           / |
//  |   2   7   3   |
//  |     \   /   8 |
//  |  5    3-------4
//  |     /   \   6 |
//  |   1   4   \   |
//  | /           \ |
//  1---------------2
string small_mesh_junction = R"CODE(
$MeshFormat
2.2 0 8
$EndMeshFormat
$Nodes
6
1 0 0 0
2 2 0 0
3 1 1 0
4 2 1 0
5 0 2 0
6 2 2 0
$EndNodes
$Elements
8
1 1 2 1 1 1 3
2 1 2 1 1 3 5
3 1 2 1 1 3 6
4 2 2 2 2 1 2 3
5 2 2 2 2 1 3 5
6 2 2 2 2 3 2 4
7 2 2 2 2 3 6 5
8 2 2 2 2 3 4 6
$EndElements
)CODE";


TEST(DOFHandler, test_all) {
  
  // distribute dofs for continuous P1 finite element.
  // The test checks that the dofs are
  // shared by adjacent elements
  // except for the dofs lying on line 1
  // without the center point X.
  {
    Mesh * mesh = mesh_constructor("{mesh_file=\"\"}", Input::FileFormat::format_JSON, MPI_COMM_WORLD);
    stringstream in(small_mesh.c_str());
    mesh->read_gmsh_from_stream(in);
    
    FE_P<1,1,3> fe1;
    FE_P<1,2,3> fe2;
    FE_P<1,3,3> fe3;
    std::shared_ptr<DiscreteSpace> ds = std::make_shared<EqualOrderDiscreteSpace>(mesh, &fe1, &fe2, &fe3);
    DOFHandlerMultiDim dh(*mesh);
    dh.distribute_dofs(ds);
    
    EXPECT_EQ( 8, dh.n_global_dofs() );
    
    unsigned int indices[5][dh.max_elem_dofs()];
    for (unsigned int i=0; i<5; i++)
      dh.get_dof_indices(mesh->element.find_id(i+1), indices[i]);
    
    // dof at node 1 is shared by elements 2, 3
    EXPECT_EQ( indices[1][0], indices[2][0] );
    
    // dof at node 2 is shared by elements 2, 4
    EXPECT_EQ( indices[1][1], indices[3][1] );
    
    // dof at node 3 is shared by elements 2, 3, 4, 5
    EXPECT_EQ( indices[4][0], indices[2][1] );
    EXPECT_EQ( indices[2][1], indices[1][2] );
    EXPECT_EQ( indices[1][2], indices[3][0] );
    
    // dof at node 3 is NOT shared by elements 1 and 5
    EXPECT_NE( indices[0][0], indices[4][0] );
    
    // dof at node 4 is shared by elements 3, 5
    EXPECT_EQ( indices[2][2], indices[4][2] );
    
    // dof at node 5 is NOT shared by elements 1, 4 and 5
    EXPECT_NE( indices[4][1], indices[0][1] );
    EXPECT_NE( indices[4][1], indices[3][2] );
    
    delete mesh;

  }

  
  
  
  
  
  // distribute dofs for continuous P1 finite element.
  // The test checks that the dofs are
  // shared by adjacent elements
  // except for the dofs separated by line elements.
  {
    Mesh * mesh = mesh_constructor("{mesh_file=\"\"}", Input::FileFormat::format_JSON, MPI_COMM_WORLD);
    stringstream in(small_mesh_junction.c_str());
    mesh->read_gmsh_from_stream(in);
    
    FE_P<1,1,3> fe1;
    FE_P<1,2,3> fe2;
    FE_P<1,3,3> fe3;
    std::shared_ptr<DiscreteSpace> ds = std::make_shared<EqualOrderDiscreteSpace>(mesh, &fe1, &fe2, &fe3);
    DOFHandlerMultiDim dh(*mesh);
    dh.distribute_dofs(ds);
    
    EXPECT_EQ( 15, dh.n_global_dofs() );
    
    unsigned int indices[mesh->n_elements()][dh.max_elem_dofs()];
    for (unsigned int i=0; i<mesh->n_elements(); i++)
      dh.get_dof_indices(mesh->element.find_id(i+1), indices[i]);
    
    // dof at node 1 is not shared by elements 1, 4, 5
    EXPECT_NE( indices[0][0], indices[3][0] );
    EXPECT_NE( indices[4][0], indices[3][0] );
    
    // dof at node 2 is shared by elements 4, 6
    EXPECT_EQ( indices[3][1], indices[5][1] );
    
    // dof at node 3 is shared by elements 4, 6, 8
    EXPECT_EQ( indices[3][2], indices[5][0] );
    EXPECT_EQ( indices[3][2], indices[7][0] );
    
    // dof at node 3 is shared by elements 1, 2, 3
    EXPECT_EQ( indices[0][1], indices[1][0] );
    EXPECT_EQ( indices[1][0], indices[2][0] );
    
    // dof at node 3 is NOT shared by elements 1, 4, 5, 7
    EXPECT_NE( indices[0][1], indices[3][2] );
    EXPECT_NE( indices[3][2], indices[4][1] );
    EXPECT_NE( indices[4][1], indices[6][0] );
    
    // dof at node 4 is shared by elements 6, 8
    EXPECT_EQ( indices[5][2], indices[7][1] );
    
    // dof at node 5 is NOT shared by elements 2, 5, 7
    EXPECT_NE( indices[1][1], indices[4][2] );
    EXPECT_NE( indices[4][2], indices[6][2] );
    
    // dof at node 6 is NOT shared by elements 3, 7, 8
    EXPECT_NE( indices[2][1], indices[6][1] );
    EXPECT_NE( indices[6][1], indices[7][2] );
    
    delete mesh;

  }

}


