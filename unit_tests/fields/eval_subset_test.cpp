/*
 * eval_subset_test.cpp
 *
 *  Created on: Dec 03, 2019
 *      Author: David Flanderka
 *
 *  Tests EvalPoints, Integral classes ...
 */

#define TEST_USE_MPI
#define FEAL_OVERRIDE_ASSERTS
#include <flow_gtest_mpi.hh>
#include <mesh_constructor.hh>
#include "arma_expect.hh"

#include "fields/eval_points.hh"
#include "fields/eval_subset.hh"
#include "fields/field_value_cache.hh"
#include "quadrature/quadrature.hh"
#include "quadrature/quadrature_lib.hh"
#include "fem/dofhandler.hh"
#include "fem/dh_cell_accessor.hh"
#include "mesh/mesh.h"
#include "system/sys_profiler.hh"


template<unsigned int dim>
arma::vec::fixed<dim> loc_coords(std::shared_ptr<EvalPoints> eval_points, unsigned int eval_point_idx) {
    return eval_points->local_point<dim>( eval_point_idx );
}


TEST(EvalPointsTest, all) {
	std::shared_ptr<EvalPoints> eval_points = std::make_shared<EvalPoints>();
	EXPECT_EQ(eval_points->size(3), 0);
	EXPECT_EQ(eval_points->n_subsets(3), 0);

    Quadrature *q_bulk_3 = new QGauss(3, 2); // dim 3
    eval_points->add_bulk<3>(*q_bulk_3 );
	EXPECT_EQ(eval_points->size(3), 4);
	EXPECT_EQ(eval_points->n_subsets(3), 1);
	EXPECT_EQ(eval_points->subset_begin(3, 0), 0);
	EXPECT_EQ(eval_points->subset_end(3, 0), 4);
	EXPECT_EQ(eval_points->subset_size(3, 0), 4);

	Quadrature *q_bulk_0 = new QGauss(0, 2); // dim 0
    eval_points->add_bulk<0>(*q_bulk_0 );
	EXPECT_EQ(eval_points->size(0), 1);
	EXPECT_EQ(eval_points->n_subsets(0), 1);
	EXPECT_EQ(eval_points->subset_begin(0, 0), 0);
	EXPECT_EQ(eval_points->subset_end(0, 0), 1);
	EXPECT_EQ(eval_points->subset_size(0, 0), 1);
}


TEST(IntegralTest, integrals_3d) {
    FilePath::set_io_dirs(".",UNIT_TESTS_SRC_DIR,"",".");
    Profiler::instance();
    PetscInitialize(0,PETSC_NULL,PETSC_NULL,PETSC_NULL);

	std::shared_ptr<EvalPoints> eval_points = std::make_shared<EvalPoints>();
    Quadrature *q_bulk = new QGauss(3, 2);
    Quadrature *q_side = new QGauss(2, 2);
    std::shared_ptr<BulkIntegral> bulk_integral = eval_points->add_bulk<3>(*q_bulk );
    std::shared_ptr<EdgeIntegral> edge_integral = eval_points->add_edge<3>(*q_side );
    std::shared_ptr<CouplingIntegral> coupling_integral = eval_points->add_coupling<3>(*q_side );
    std::shared_ptr<BoundaryIntegral> boundary_integral = eval_points->add_boundary<3>(*q_side );

    Mesh * mesh = mesh_full_constructor("{ mesh_file=\"mesh/simplest_cube.msh\", optimize_mesh=false }");
    std::shared_ptr<DOFHandlerMultiDim> dh = std::make_shared<DOFHandlerMultiDim>(*mesh);
    DHCellAccessor dh_cell(dh.get(), 3);
    ElementCacheMap elm_cache_map;
    elm_cache_map.init(eval_points);

    {
        // Test of bulk local points
    	std::vector<arma::vec3> expected_vals = {{0.138196601125010504, 0.138196601125010504, 0.138196601125010504},
    			                                 {0.138196601125010504, 0.138196601125010504, 0.585410196624968515},
												 {0.138196601125010504, 0.585410196624968515, 0.138196601125010504},
												 {0.585410196624968515, 0.138196601125010504, 0.138196601125010504}};
    	unsigned int i=0; // iter trought expected_vals
    	for (auto p : bulk_integral->points(dh_cell, &elm_cache_map)) {
            EXPECT_ARMA_EQ(loc_coords<3>(eval_points, p.eval_point_idx()), expected_vals[i]);
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
        expected_vals[2] = { {0.0, 0.666666666666666657, 0.166666666666666657},
	                         {0.0, 0.166666666666666657, 0.166666666666666741},
                             {0.0, 0.166666666666666741, 0.666666666666666657} };
        expected_vals[3] = { {0.666666666666666741, 0.166666666666666657, 0.166666666666666657},
                             {0.166666666666666657, 0.166666666666666657, 0.666666666666666741},
                             {0.166666666666666657, 0.666666666666666741, 0.166666666666666657} };
        unsigned int i_side=0, i_point; // iter trought expected_vals
        for (auto side_acc : dh_cell.side_range()) {
            i_point=0;
            for ( auto p : edge_integral->points(side_acc, &elm_cache_map) ) {
            	EXPECT_ARMA_EQ(loc_coords<3>(eval_points, p.eval_point_idx()), expected_vals[i_side][i_point]);
                ++i_point;
            }
            ++i_side;
        }
    }

    {
        // Test of boundary
        std::vector< std::vector<arma::vec> > expected_vals(3);
        expected_vals[0] = { {0.0, 0.666666666666666657, 0.166666666666666657},
	                         {0.0, 0.166666666666666657, 0.166666666666666741},
                             {0.0, 0.166666666666666741, 0.666666666666666657} };
        expected_vals[1] = { {0.666666666666666741, 0.166666666666666657, 0.166666666666666657},
                             {0.166666666666666657, 0.166666666666666657, 0.666666666666666741},
                             {0.166666666666666657, 0.666666666666666741, 0.166666666666666657} };
        expected_vals[2] = { {0.166666666666666657, 0.166666666666666657},
                             {0.166666666666666657, 0.666666666666666741},
                             {0.666666666666666741, 0.166666666666666657} };
        unsigned int i_side=0, i_point; // iter trought expected_vals
        for (auto side_acc : dh_cell.side_range()) {
            if (! side_acc.side().is_boundary())
                continue;
            i_point=0;
            for ( auto p : boundary_integral->points(side_acc, &elm_cache_map) ) {
                EXPECT_ARMA_EQ(loc_coords<3>(eval_points, p.eval_point_idx()), expected_vals[i_side][i_point]);
                auto p_bdr = p.point_bdr( side_acc.cond().element_accessor() );
                EXPECT_ARMA_EQ(loc_coords<2>(eval_points, p_bdr.eval_point_idx()), expected_vals[2][i_point]);
                ++i_point;
            }
            ++i_side;
        }
    }

    {
        // Test of neighbours (3D - 2D - 3D elements)
        DHCellAccessor dh_ngh_cell(dh.get(), 1);
        std::vector< std::vector<arma::vec> > expected_vals(3);
        expected_vals[0] = { {0.166666666666666657, 0.166666666666666657},         // local points on cell (element of lower dim)
                             {0.166666666666666657, 0.666666666666666741},
                             {0.666666666666666741, 0.166666666666666657} };
        expected_vals[1] = { {0.166666666666666657, 0.0, 0.666666666666666741},    // local points on side (first element of higher dim)
                             {0.666666666666666741, 0.0, 0.166666666666666657},
                             {0.166666666666666657, 0.0, 0.166666666666666657} };
        expected_vals[2] = { {0.166666666666666657, 0.666666666666666741, 0.0},    // local points on side (second element of higher dim)
                             {0.666666666666666741, 0.166666666666666657, 0.0},
                             {0.166666666666666657, 0.166666666666666657, 0.0} };
        unsigned int i_side=1, i_point; // iterates trought expected_vals
        for (auto ngh_side_acc : dh_ngh_cell.neighb_sides()) {
            i_point=0;
            for ( auto p_side : coupling_integral->points(ngh_side_acc, &elm_cache_map) ) {
            	auto p_cell = p_side.lower_dim(dh_ngh_cell);
                EXPECT_ARMA_EQ(loc_coords<2>(eval_points, p_cell.eval_point_idx()), expected_vals[0][i_point]);
                EXPECT_ARMA_EQ(loc_coords<3>(eval_points, p_side.eval_point_idx()), expected_vals[i_side][i_point]);
                ++i_point;
            }
            ++i_side;
        }
    }
}

