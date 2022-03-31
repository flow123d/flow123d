#define TEST_USE_PETSC
#include <flow_gtest_mpi.hh>
#include <cmath>
#include "mesh/bc_mesh.hh"
#include "mesh/accessors.hh"
#include "mesh/neighbours.h"
#include <mesh_constructor.hh>



void print_bc_mesh(BCMesh *bcmesh)
{
  for (auto el : bcmesh->elements_range())
  {
    cout << "element " << el.idx() << " dim " << el->dim() << " nodes";
    for (unsigned int nid=0; nid<el->n_nodes(); nid++)
      cout << " " << el.node(nid).index();
    cout << endl;
  }

  for (auto edge : bcmesh->edge_range())
  {
    cout << "edge " << edge.idx() << " elements";
    for (unsigned int sid=0; sid<edge.n_sides(); sid++)
      cout << " " << edge.side(sid)->elem_idx();
    cout << endl;
  }

  for (auto nbid = 0; nbid<bcmesh->n_vb_neighbours(); nbid++)
  {
    auto nb = bcmesh->vb_neighbour(nbid);
    cout << "neighbour " << nbid << " edge " << nb.edge_idx() << " elements " << nb.side()->elem_idx() << " " << nb.element().idx() << endl;
  }
}


TEST(BCMesh, test_all) {
  
  // Create mesh for boundary elements.
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
    Mesh * mesh = mesh_full_constructor("{ mesh_file=\"fem/small_mesh.msh\", optimize_mesh=false }");
    BCMesh * bcmesh = mesh->bc_mesh();

    print_bc_mesh(bcmesh);

    EXPECT_EQ( 6, bcmesh->n_elements() );
    EXPECT_EQ( 5, bcmesh->n_edges() );
    EXPECT_EQ( 2, bcmesh->n_vb_neighbours() );

    // check element dimensions
    vector<int> elem_dims(4, 0);
    for (auto el : bcmesh->elements_range()) elem_dims[el->dim()]++;
    EXPECT_EQ( 2, elem_dims[0] );
    EXPECT_EQ( 4, elem_dims[1] );

    // check neighbours
    for (unsigned int nbid=0; nbid<bcmesh->n_vb_neighbours(); nbid++)
    {
      EXPECT_EQ( 0, bcmesh->vb_neighbour(nbid).element()->dim() );
      EXPECT_EQ( 1, bcmesh->vb_neighbour(nbid).edge().n_sides() );
    }
    
    delete mesh;
  }

  
  
  
  
  // Create mesh for boundary elements.
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
    Mesh * mesh = mesh_full_constructor("{ mesh_file=\"fem/small_mesh_junction.msh\", optimize_mesh=false }");
    BCMesh * bcmesh = mesh->bc_mesh();

    print_bc_mesh(bcmesh);

    EXPECT_EQ( 8, bcmesh->n_elements() );
    EXPECT_EQ( 8, bcmesh->n_edges() );
    EXPECT_EQ( 6, bcmesh->n_vb_neighbours() );

    // check element dimensions
    vector<int> elem_dims(4, 0);
    for (auto el : bcmesh->elements_range()) elem_dims[el->dim()]++;
    EXPECT_EQ( 3, elem_dims[0] );
    EXPECT_EQ( 5, elem_dims[1] );

    // check neighbours
    for (unsigned int nbid=0; nbid<bcmesh->n_vb_neighbours(); nbid++)
    {
      EXPECT_EQ( 0, bcmesh->vb_neighbour(nbid).element()->dim() );
      EXPECT_EQ( 1, bcmesh->vb_neighbour(nbid).edge().n_sides() );
    }
    
    delete mesh;

  }
}


