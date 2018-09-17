#define TEST_USE_PETSC
#include <flow_gtest_mpi.hh>
#include <cmath>
#include "fem/fe_p.hh"
#include "mesh/mesh.h"
#include <mesh_constructor.hh>
#include "fem/dofhandler.hh"





TEST(DOFHandler, test_all) {
  
  // distribute dofs for continuous P1 finite element.
  // The test checks that the dofs are
  // shared by adjacent elements
  // except for the dofs lying on line 1
  // without the center point 3.
  {
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
    FilePath::set_io_dirs(".",UNIT_TESTS_SRC_DIR,"",".");
    Mesh * mesh = mesh_full_constructor("{mesh_file=\"fem/small_mesh.msh\"}");
    
    FE_P<0> fe0(1);
    FE_P<1> fe1(1);
    FE_P<2> fe2(1);
    FE_P<3> fe3(1);
    std::shared_ptr<DiscreteSpace> ds = std::make_shared<EqualOrderDiscreteSpace>(mesh, &fe0, &fe1, &fe2, &fe3);
    DOFHandlerMultiDim dh(*mesh);
    dh.distribute_dofs(ds, true);
    
    EXPECT_EQ( 8, dh.n_global_dofs() );
    
    std::vector<int> indices[5];
    for (unsigned int i=0; i<5; i++)
    {
      indices[i].resize(dh.max_elem_dofs());
      dh.get_dof_indices(mesh->element_accessor(i), indices[i]);
    }
    
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
    // simple mesh consisting of 5 triangles and 3 line elements
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
    FilePath::set_io_dirs(".",UNIT_TESTS_SRC_DIR,"",".");
    Mesh * mesh = mesh_full_constructor("{mesh_file=\"fem/small_mesh_junction.msh\"}");
    
    FE_P<0> fe0(1);
    FE_P<1> fe1(1);
    FE_P<2> fe2(1);
    FE_P<3> fe3(1);
    std::shared_ptr<DiscreteSpace> ds = std::make_shared<EqualOrderDiscreteSpace>(mesh, &fe0, &fe1, &fe2, &fe3);
    DOFHandlerMultiDim dh(*mesh);
    dh.distribute_dofs(ds, true);
    
    EXPECT_EQ( 15, dh.n_global_dofs() );
    
    std::vector<int> indices[mesh->n_elements()];
    for (unsigned int i=0; i<mesh->n_elements(); i++)
    {
      indices[i].resize(dh.max_elem_dofs());
      dh.get_dof_indices(mesh->element_accessor(i), indices[i]);
    }
    
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


