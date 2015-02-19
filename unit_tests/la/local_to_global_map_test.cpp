/*
 * local_to_global_map.cc
 *
 *  Created on: Mar 9, 2012
 *      Author: jb
 */



#define TEST_USE_PETSC

#include <flow_gtest_mpi.hh>

#include <la/local_to_global_map.hh>
#include <la/distribution.hh>




TEST(la, local_to_global_map) {

    // local setup
    int np, rank;
    MPI_Comm_rank(PETSC_COMM_WORLD, &rank);
    MPI_Comm_size(PETSC_COMM_WORLD, &np);

    if (np > 1)
    {
    Distribution my_distr(3+rank, PETSC_COMM_WORLD);
    LocalToGlobalMap lg_map(my_distr);

    lg_map.insert( my_distr.begin() );
    lg_map.insert( my_distr.begin() + 2 );
    lg_map.insert( (my_distr.end() +1) % my_distr.size() );
    lg_map.insert( (my_distr.end() +1) % my_distr.size() );
    EXPECT_EQ( 0, lg_map.size() );

    lg_map.finalize();

    // test size
    EXPECT_EQ(my_distr.lsize() + 1,  lg_map.size() );
    // test operator[]
    EXPECT_EQ( my_distr.begin(),       lg_map[0]);  // first local

    EXPECT_EQ( my_distr.end() - 1 - rank,      lg_map[2]); // last local

    EXPECT_EQ( (my_distr.end() +1) % my_distr.size() ,      lg_map[my_distr.lsize()]); // nonlocal

    // test inner vector
    EXPECT_EQ(my_distr.lsize() + 1,  lg_map.get_map_vector().size() );

    EXPECT_EQ( my_distr.begin(),       lg_map.get_map_vector()[0]);

    EXPECT_EQ( my_distr.end() - 1 - rank,      lg_map.get_map_vector()[2]);

    EXPECT_EQ((my_distr.end() +1) % my_distr.size(), lg_map.get_map_vector()[my_distr.lsize()]);

    // test distribution getter
    boost::shared_ptr<Distribution> sh_distr=boost::make_shared<Distribution>(3+rank,PETSC_COMM_WORLD);
    LocalToGlobalMap sh_lg_map(sh_distr);
    EXPECT_TRUE( sh_distr == sh_lg_map.get_distr() );
    EXPECT_EQ( 2, sh_distr.use_count() );
    EXPECT_EQ( 2, sh_lg_map.get_distr().use_count() );

    }
}


