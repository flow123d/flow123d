/*
 * eval_subset_test.cpp
 *
 *  Created on: Dec 03, 2019
 *      Author: David Flanderka
 *
 *  Tests EvalPoints, EvalSubset, BulkPoint and SidePoint classes
 */

#define TEST_USE_MPI
#define FEAL_OVERRIDE_ASSERTS
#include <flow_gtest_mpi.hh>
#include <mesh_constructor.hh>
#include "arma_expect.hh"

#include "fields/eval_points.hh"
#include "fields/eval_subset.hh"
#include "quadrature/quadrature.hh"
#include "quadrature/quadrature_lib.hh"
#include "fem/dofhandler.hh"
#include "fem/dh_cell_accessor.hh"
#include "mesh/mesh.h"
#include "mesh/sides.h"
#include "mesh/side_impl.hh"
#include "system/sys_profiler.hh"


TEST(EvalPointsTest, all) {
	std::shared_ptr<EvalPoints> eval_points = std::make_shared<EvalPoints>();
	EXPECT_EQ(eval_points->point_dim(), EvalPoints::undefined_dim);
	EXPECT_EQ(eval_points->size(), 0);
	EXPECT_EQ(eval_points->n_subsets(), 0);

    Quadrature *q_bulk = new QGauss(3, 2);
    eval_points->add_bulk<3>(*q_bulk );
	EXPECT_EQ(eval_points->point_dim(), 3);
	EXPECT_EQ(eval_points->size(), 4);
	EXPECT_EQ(eval_points->n_subsets(), 1);
	EXPECT_EQ(eval_points->subset_begin(0), 0);
	EXPECT_EQ(eval_points->subset_end(0), 4);
	EXPECT_EQ(eval_points->subset_size(0), 4);
}


TEST(EvalSubsetTest, subsets_3d) {
    FilePath::set_io_dirs(".",UNIT_TESTS_SRC_DIR,"",".");
    Profiler::initialize();
    PetscInitialize(0,PETSC_NULL,PETSC_NULL,PETSC_NULL);

	std::shared_ptr<EvalPoints> eval_points = std::make_shared<EvalPoints>();
    Quadrature *q_bulk = new QGauss(3, 2);
    Quadrature *q_side = new QGauss(2, 2);
    std::shared_ptr<EvalSubset> bulk_points = eval_points->add_bulk<3>(*q_bulk );
    std::shared_ptr<EvalSubset> side_points = eval_points->add_side<3>(*q_side );


    Mesh * mesh = mesh_full_constructor("{mesh_file=\"mesh/simplest_cube.msh\"}");
    std::shared_ptr<DOFHandlerMultiDim> dh = std::make_shared<DOFHandlerMultiDim>(*mesh);
    DHCellAccessor dh_cell(dh.get(), 3);

    {
        // Test of bulk local points
    	std::vector<arma::vec3> expected_vals = {{0.138196601125010504, 0.138196601125010504, 0.138196601125010504},
    			                                 {0.138196601125010504, 0.138196601125010504, 0.585410196624968515},
												 {0.138196601125010504, 0.585410196624968515, 0.138196601125010504},
												 {0.585410196624968515, 0.138196601125010504, 0.138196601125010504}};
    	unsigned int i=0; // iter trought expected_vals
    	for (auto p : bulk_points->points(dh_cell)) {
            EXPECT_ARMA_EQ(p.loc_coords<3>(), expected_vals[i]);
			++i;
        }
    }
    {
        // Test of side local points
        std::vector< std::vector<arma::vec3> > expected_vals(4);
        expected_vals[0] = { {0.166666666666666657, 0.166666666666666657, 0.0},
                             {0.666666666666666741, 0.166666666666666657, 0.0},
                             {0.166666666666666657, 0.666666666666666741, 0.0} };
        expected_vals[1] = { {0.166666666666666657, 0.0, 0.166666666666666657},
                             {0.166666666666666657, 0.0, 0.666666666666666741},
                             {0.666666666666666741, 0.0, 0.166666666666666657} };
        expected_vals[2] = { {0.0, 0.166666666666666657, 0.166666666666666657},
	                         {0.0, 0.166666666666666657, 0.666666666666666741},
                             {0.0, 0.666666666666666741, 0.166666666666666657} };
        expected_vals[3] = { {0.666666666666666741, 0.166666666666666657, 0.166666666666666657},
                             {0.166666666666666657, 0.166666666666666657, 0.666666666666666741},
                             {0.166666666666666657, 0.666666666666666741, 0.166666666666666657} };
        unsigned int i_side=0, i_point; // iter trought expected_vals
        for (auto side_acc : dh_cell.side_range()) {
        	i_point=0;
            for ( auto p : side_points->points(side_acc) ) {
            	EXPECT_ARMA_EQ(p.loc_coords<3>(), expected_vals[i_side][i_point]);
                ++i_point;
            }
            ++i_side;
        }
    }
}

