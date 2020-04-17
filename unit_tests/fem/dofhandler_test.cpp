#define TEST_USE_PETSC
#include <flow_gtest_mpi.hh>
#include <cmath>
#include "fem/fe_p.hh"
#include "fem/fe_rt.hh"
#include "fem/fe_system.hh"
#include "mesh/mesh.h"
#include <mesh_constructor.hh>
#include "fem/dofhandler.hh"
#include "fem/dh_cell_accessor.hh"
#include "tools/mixed.hh"





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
    
    MixedPtr<FE_P> fe(1);
    std::shared_ptr<DiscreteSpace> ds = std::make_shared<EqualOrderDiscreteSpace>(mesh, fe);
    DOFHandlerMultiDim dh(*mesh);
    dh.distribute_dofs(ds);
    
    EXPECT_EQ( 8, dh.n_global_dofs() );
    
    dh.print();
    
    std::vector<int> indices[5];
    std::vector<bool> own_elem(5, false); // hold if cell is own or ghost on process and can be evaluated (see tests)
    for ( DHCellAccessor cell : dh.local_range() )
    {
        auto elem_idx = cell.elm_idx();
        own_elem[elem_idx] = true;
        indices[elem_idx].resize(dh.max_elem_dofs());
    	cell.get_dof_indices(indices[elem_idx]);
    }

    // dof at node 1 is shared by elements 2, 3
    if (own_elem[1] & own_elem[2]) EXPECT_EQ( indices[1][0], indices[2][0] );
    
    // dof at node 2 is shared by elements 2, 4
    if (own_elem[1] & own_elem[3]) EXPECT_EQ( indices[1][1], indices[3][1] );
    
    // dof at node 3 is shared by elements 2, 3, 4, 5
    if (own_elem[4] & own_elem[2]) EXPECT_EQ( indices[4][0], indices[2][1] );
    if (own_elem[2] & own_elem[1]) EXPECT_EQ( indices[2][1], indices[1][2] );
    if (own_elem[1] & own_elem[3]) EXPECT_EQ( indices[1][2], indices[3][0] );
    
    // dof at node 3 is NOT shared by elements 1 and 5
    if (own_elem[0] & own_elem[4]) EXPECT_NE( indices[0][0], indices[4][0] );
    
    // dof at node 4 is shared by elements 3, 5
    if (own_elem[2] & own_elem[4]) EXPECT_EQ( indices[2][2], indices[4][2] );
    
    // dof at node 5 is NOT shared by elements 1, 4 and 5
    if (own_elem[4] & own_elem[0]) EXPECT_NE( indices[4][1], indices[0][1] );
    if (own_elem[4] & own_elem[3]) EXPECT_NE( indices[4][1], indices[3][2] );
    
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
    
    MixedPtr<FE_P> fe(1);
    std::shared_ptr<DiscreteSpace> ds = std::make_shared<EqualOrderDiscreteSpace>(mesh, fe);
    DOFHandlerMultiDim dh(*mesh);
    dh.distribute_dofs(ds);
    
    EXPECT_EQ( 15, dh.n_global_dofs() );
    
    dh.print();
    
    std::vector<int> indices[8];
    std::vector<bool> own_elem(8, false); // hold if cell is own or ghost on process and can be evaluated (see tests)
    for ( DHCellAccessor cell : dh.local_range() )
    {
        auto elem_idx = cell.elm_idx();
        own_elem[elem_idx] = true;
        indices[elem_idx].resize(dh.max_elem_dofs());
        cell.get_dof_indices(indices[elem_idx]);
    }
    
    // dof at node 1 is not shared by elements 1, 4, 5
    if (own_elem[0] & own_elem[3]) EXPECT_NE( indices[0][0], indices[3][0] );
    if (own_elem[4] & own_elem[3]) EXPECT_NE( indices[4][0], indices[3][0] );
    
    // dof at node 2 is shared by elements 4, 6
    if (own_elem[3] & own_elem[5]) EXPECT_EQ( indices[3][1], indices[5][1] );
    
    // dof at node 3 is shared by elements 4, 6, 8
    if (own_elem[3] & own_elem[5]) EXPECT_EQ( indices[3][2], indices[5][0] );
    if (own_elem[3] & own_elem[7]) EXPECT_EQ( indices[3][2], indices[7][0] );
    
    // dof at node 3 is shared by elements 1, 2, 3
    if (own_elem[0] & own_elem[1]) EXPECT_EQ( indices[0][1], indices[1][0] );
    if (own_elem[1] & own_elem[2]) EXPECT_EQ( indices[1][0], indices[2][0] );
    
    // dof at node 3 is NOT shared by elements 1, 4, 5, 7
    if (own_elem[0] & own_elem[3]) EXPECT_NE( indices[0][1], indices[3][2] );
    if (own_elem[3] & own_elem[4]) EXPECT_NE( indices[3][2], indices[4][1] );
    if (own_elem[4] & own_elem[6]) EXPECT_NE( indices[4][1], indices[6][0] );
    
    // dof at node 4 is shared by elements 6, 8
    if (own_elem[5] & own_elem[7]) EXPECT_EQ( indices[5][2], indices[7][1] );
    
    // dof at node 5 is NOT shared by elements 2, 5, 7
    if (own_elem[1] & own_elem[4]) EXPECT_NE( indices[1][1], indices[4][2] );
    if (own_elem[4] & own_elem[6]) EXPECT_NE( indices[4][2], indices[6][2] );
    
    // dof at node 6 is NOT shared by elements 3, 7, 8
    if (own_elem[2] & own_elem[6]) EXPECT_NE( indices[2][1], indices[6][1] );
    if (own_elem[6] & own_elem[7]) EXPECT_NE( indices[6][1], indices[7][2] );
    
    delete mesh;

  }

}


TEST(DOFHandler, test_sub_handler)
{
    FESystem<0> fe_sys0({ std::make_shared<FE_P_disc<0> >(0),
                          std::make_shared<FE_P<0> >(0),
                          std::make_shared<FE_CR<0> >() });
    FESystem<1> fe_sys1({ std::make_shared<FE_RT0<1> >(),
                          std::make_shared<FE_P<1> >(0),
                          std::make_shared<FE_CR<1> >() });
    FESystem<2> fe_sys2({ std::make_shared<FE_RT0<2> >(),
                          std::make_shared<FE_P<2> >(0),
                          std::make_shared<FE_CR<2> >() });
    FESystem<3> fe_sys3({ std::make_shared<FE_RT0<3> >(),
                          std::make_shared<FE_P<3> >(0),
                          std::make_shared<FE_CR<3> >() });
    MixedPtr<FESystem> fe_sys( std::make_shared<FESystem<0>>(fe_sys0), std::make_shared<FESystem<1>>(fe_sys1),
    	                                    std::make_shared<FESystem<2>>(fe_sys2), std::make_shared<FESystem<3>>(fe_sys3) );

    FilePath::set_io_dirs(".",UNIT_TESTS_SRC_DIR,"",".");
    Mesh * mesh = mesh_full_constructor("{mesh_file=\"fem/small_mesh_junction.msh\"}");
    std::shared_ptr<DiscreteSpace> ds = std::make_shared<EqualOrderDiscreteSpace>(mesh, fe_sys);
    std::shared_ptr<DOFHandlerMultiDim> dh = std::make_shared<DOFHandlerMultiDim>(*mesh);
    dh->distribute_dofs(ds);
    std::shared_ptr<SubDOFHandlerMultiDim> sub_dh = std::make_shared<SubDOFHandlerMultiDim>(dh, 0);
    std::vector<LocDofVec > loc_indices(mesh->n_elements());
    LocDofVec loc_sub_indices;

    dh->print();    
    sub_dh->print();
    
    VectorMPI vec = dh->create_vector();
    VectorMPI subvec = sub_dh->create_vector();
    
    // init cell dof indices
    for (auto cell : dh->local_range())
        loc_indices[cell.elm_idx()] = cell.get_loc_dof_indices();
    
    // init vec and update subvec
    for (auto cell : dh->own_range())
        for (unsigned int i=0; i<dh->ds()->n_elem_dofs(cell.elm()); i++)
            vec[loc_indices[cell.elm_idx()][i]] = cell.elm_idx()*dh->max_elem_dofs()+i;
    vec.local_to_ghost_begin();
    vec.local_to_ghost_end();
    sub_dh->update_subvector(vec, subvec);
    
    // check that dofs on sub_dh are equal to dofs on dh
    for (auto cell : sub_dh->local_range())
    {
        loc_sub_indices = cell.get_loc_dof_indices();
        for (unsigned int i=0; i<cell.n_dofs(); i++)
        {
            // local indices
            EXPECT_EQ( sub_dh->parent_indices()[loc_sub_indices[i]], loc_indices[cell.elm_idx()][i] );
            // values in mpi vectors
            EXPECT_EQ( vec[loc_indices[cell.elm_idx()][i]], subvec[loc_sub_indices[i]] );
        }
    }

    // modify subvec and update "parent" vec
    for (auto cell : sub_dh->own_range())
    {
        loc_sub_indices = cell.get_loc_dof_indices();
        for (unsigned int i=0; i<sub_dh->ds()->n_elem_dofs(cell.elm()); i++)
            subvec[loc_sub_indices[i]] = -(cell.elm_idx()*dh->max_elem_dofs()+i);
    }
    subvec.local_to_ghost_begin();
    subvec.local_to_ghost_end();
    sub_dh->update_parent_vector(vec, subvec);
    // check values in mpi vectors
    for (auto cell : sub_dh->local_range())
    {
        loc_sub_indices = cell.get_loc_dof_indices();
        for (unsigned int i=0; i<cell.n_dofs(); i++)
            EXPECT_EQ( vec[loc_indices[cell.elm_idx()][i]], subvec[loc_sub_indices[i]] );
    }
    
    delete mesh;
    
}



// distribute dofs for continuous RT0 finite element.
// The test checks that the dofs are
// shared by adjacent elements
// except for the dofs separated by line elements.
TEST(DOFHandler, test_rt)
{
    FilePath::set_io_dirs(".",UNIT_TESTS_SRC_DIR,"",".");
    Mesh * mesh = mesh_full_constructor("{mesh_file=\"fem/small_mesh_junction.msh\"}");

    MixedPtr<FE_RT0> fe;
    std::shared_ptr<DiscreteSpace> ds = std::make_shared<EqualOrderDiscreteSpace>(mesh, fe);
    DOFHandlerMultiDim dh(*mesh);
    dh.distribute_dofs(ds);

     EXPECT_EQ( 17, dh.n_global_dofs() );

    dh.print();

    std::vector<int> indices[8];
    std::vector<bool> own_elem(8, false); // hold if cell is own or ghost on process and can be evaluated (see tests)
    for ( DHCellAccessor cell : dh.local_range() )
    {
        auto elem_idx = cell.elm_idx();
        own_elem[elem_idx] = true;
        indices[elem_idx].resize(dh.max_elem_dofs());
        cell.get_dof_indices(indices[elem_idx]);
    }

    // dof at el. 1 side 1 equals dof at el. 2 side 0
    if (own_elem[0] & own_elem[1]) EXPECT_EQ( indices[0][1], indices[1][0] );
    
    // dof at el. 1 side 1 equals dof at el. 3 side 0
    if (own_elem[0] & own_elem[2]) EXPECT_EQ( indices[0][1], indices[2][0] );
    
    // dof at el. 4 side 2 equals dof at el. 6 side 0
    if (own_elem[3] & own_elem[5]) EXPECT_EQ( indices[3][2], indices[5][0] );
    
    // dof at el. 6 side 1 equals dof at el. 8 side 0
    if (own_elem[5] & own_elem[7]) EXPECT_EQ( indices[5][1], indices[7][0] );

    delete mesh;

}



TEST(DHAccessors, dh_cell_accessors) {
    FilePath::set_io_dirs(".",UNIT_TESTS_SRC_DIR,"",".");
    Mesh * mesh = mesh_full_constructor("{mesh_file=\"mesh/simplest_cube.msh\"}");

    DOFHandlerMultiDim dh(*mesh);
    {
        MixedPtr<FE_P> fe(1);
        std::shared_ptr<DiscreteSpace> ds = std::make_shared<EqualOrderDiscreteSpace>(mesh, fe);
        dh.distribute_dofs(ds);
    }
    DOFHandlerMultiDim dh_2(*mesh); // test cell_with_other_dh method
    {
        MixedPtr<FE_RT0> fe;
        std::shared_ptr<DiscreteSpace> ds = std::make_shared<EqualOrderDiscreteSpace>(mesh, fe);
        dh_2.distribute_dofs(ds);
    }
    unsigned int i_distr=0;

    std::vector<unsigned int> side_elm_idx, neigh_elem_idx;
    for( DHCellAccessor cell : dh.own_range() ) {
    	EXPECT_EQ( cell.elm_idx(), dh.mesh()->get_el_4_loc()[i_distr] );
    	DHCellAccessor cell_2 = cell.cell_with_other_dh(&dh_2);
    	EXPECT_EQ( cell.local_idx(), cell_2.local_idx() );

    	for( DHCellSide cell_side : cell.side_range() ) {
            EXPECT_EQ( cell.elm_idx(), cell_side.elem_idx() );
        	side_elm_idx.clear();
        	for( DHCellSide edge_side : cell_side.edge_sides() ) {
        		side_elm_idx.push_back( edge_side.elem_idx() );
        	}
            Edge edg = cell_side.side().edge();
            EXPECT_EQ( side_elm_idx.size(), edg.n_sides());
            for (uint sid=0; sid<edg.n_sides(); sid++) {
            	EXPECT_EQ( side_elm_idx[sid], edg.side(sid)->element().idx());
            }
        }

        neigh_elem_idx.clear();
        for( DHCellSide neighb_side : cell.neighb_sides() ) {
        	neigh_elem_idx.push_back( neighb_side.elem_idx() );
        }
        EXPECT_EQ( neigh_elem_idx.size(), cell.elm()->n_neighs_vb());
        for (uint nid=0; nid<cell.elm()->n_neighs_vb(); nid++) {
        	EXPECT_EQ( neigh_elem_idx[nid], cell.elm()->neigh_vb[nid]->side()->elem_idx() );
        }

    	++i_distr;
    }

    delete mesh;
}
