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
    	std::vector<arma::vec3> expected_vals = {
          {0.138196601125010504, 0.138196601125010504, 0.138196601125010504},
          {0.138196601125010504, 0.138196601125010504, 0.585410196624968515},
          {0.138196601125010504, 0.585410196624968515, 0.138196601125010504},
          {0.585410196624968515, 0.138196601125010504, 0.138196601125010504}};
          
    	unsigned int i=0; // iter trought expected_vals
    	for (auto p : bulk_integral->points(elm_cache_map.position_in_cache(dh_cell.elm_idx()), &elm_cache_map)) {
    	    DebugOut() << "BULK i: " << i << "eval_point_idx: " << p.eval_point_idx();
            EXPECT_ARMA_EQ(expected_vals[i], loc_coords<3>(eval_points, p.eval_point_idx()));
			++i;
        }
    }
    {
        // Test of side local points
        // results corresponds to reference element, permutation zero
        std::vector< std::vector<arma::vec3> > expected_vals(4);
        expected_vals[0] = { {1.0/6, 1.0/6, 0.0},
                             {1.0/6, 2.0/3, 0.0},
                             {2.0/3, 1.0/6, 0.0}};
        expected_vals[1] = { {1.0/6, 0.0, 1.0/6},
                             {1.0/6, 0.0, 2.0/3},
                             {2.0/3, 0.0, 1.0/6} };
        expected_vals[2] = { {0.0, 1.0/6, 1.0/6},
                             {0.0, 1.0/6, 2.0/3},
                             {0.0, 2.0/3, 1.0/6}};
        expected_vals[3] = { {2.0/3, 1.0/6, 1.0/6},
                             {1.0/6, 1.0/6, 2.0/3},
                             {1.0/6, 2.0/3, 1.0/6} };
        unsigned int i_side=0, i_point; // iter trought expected_vals
        for (auto side_acc : dh_cell.side_range()) {
            i_point=0;
            for ( EdgePoint p : edge_integral->points(side_acc, &elm_cache_map) ) {
                DebugOut() << "side: " << i_side << " ip: " << i_point << "eval point idx: " << p.eval_point_idx();
            	EXPECT_ARMA_EQ(expected_vals[i_side][i_point], loc_coords<3>(eval_points, p.eval_point_idx()));
                ++i_point;
            }
            ++i_side;
        }
    }

    {
        // Test of boundary
        std::vector< std::vector<arma::vec> > expected_vals(3);
        expected_vals[0] = { {1.0/6, 1.0/6, 0.0},
                             {1.0/6, 2.0/3, 0.0},
                             {2.0/3, 1.0/6, 0.0}};
        expected_vals[1] = { {1.0/6, 0.0, 1.0/6},
                             {1.0/6, 0.0, 2.0/3},
                             {2.0/3, 0.0, 1.0/6} };

        expected_vals[2] = { {1.0/6, 1.0/6},
                             {1.0/6, 2.0/3},
                             {2.0/3, 1.0/6} };
        unsigned int i_side=0, i_point; // iter trought expected_vals
        for (auto side_acc : dh_cell.side_range()) {
            if (! side_acc.side().is_boundary())
                continue;
            i_point=0;
            for ( auto p : boundary_integral->points(side_acc, &elm_cache_map) ) {
                DebugOut() << "side: " << i_side << " ip: " << i_point;
                EXPECT_ARMA_EQ(expected_vals[i_side][i_point], loc_coords<3>(eval_points, p.eval_point_idx()));


                auto p_bdr = p.point_bdr( side_acc.cond().element_accessor() );
                DebugOut() << "side: " << i_side << " ip: " << i_point << "eval point idx: " << p_bdr.eval_point_idx();
                EXPECT_ARMA_EQ(expected_vals[2][i_point], loc_coords<2>(eval_points, p_bdr.eval_point_idx()));
                ++i_point;
            }
            ++i_side;
        }
    }

    {
        // Test of neighbours (3D - 2D - 3D elements)
        DHCellAccessor dh_ngh_cell(dh.get(), 1);
        std::vector< std::vector<arma::vec> > expected_vals(3);
        expected_vals[0] = { {1.0/6, 1.0/6},         // local points on cell (element of lower dim)
                             {1.0/6, 2.0/3},
                             {2.0/3, 1.0/6} };
        expected_vals[1] = { {1.0/6, 1.0/6, 0.0},   // local points on side 0 of a 3d element
                             {1.0/6, 2.0/3, 0.0},
                             {2.0/3, 1.0/6, 0.0}};
        expected_vals[2] = { {0.0, 1.0/6, 1.0/6},   // local points on side 2 of a 3d element
                             {0.0, 1.0/6, 2.0/3},
                             {0.0, 2.0/3, 1.0/6}};

        unsigned int i_side=1, i_point; // iterates trought expected_vals
        for (auto ngh_side_acc : dh_ngh_cell.neighb_sides()) {
            i_point=0;
            for ( auto p_side : coupling_integral->points(ngh_side_acc, &elm_cache_map) ) {
            	auto p_cell = p_side.lower_dim(dh_ngh_cell);
            	DebugOut() << "ele side: " << ngh_side_acc.side_idx() << " ip: " << i_point;
                EXPECT_ARMA_EQ(expected_vals[0][i_point], loc_coords<2>(eval_points, p_cell.eval_point_idx()));

                EXPECT_ARMA_EQ(expected_vals[i_side][i_point], loc_coords<3>(eval_points, p_side.eval_point_idx()));
                ++i_point;
            }
            ++i_side;
        }
    }
}


#ifdef FLOW123D_RUN_UNIT_BENCHMARKS

static const unsigned int profiler_loops = 1e7;

class RangesSpeedTest : public testing::Test, public ElementCacheMap {
public:
	RangesSpeedTest() {
        FilePath::set_io_dirs(".",UNIT_TESTS_SRC_DIR,"",".");
        Profiler::instance();
        PetscInitialize(0,PETSC_NULL,PETSC_NULL,PETSC_NULL);

    	eval_points_ = std::make_shared<EvalPoints>();
        q_bulk_ = new QGauss(3, 2);
        q_side_ = new QGauss(2, 2);
        bulk_integral_ = eval_points_->add_bulk<3>(*q_bulk_ );
        edge_integral_ = eval_points_->add_edge<3>(*q_side_ );
        coupling_integral_ = eval_points_->add_coupling<3>(*q_side_ );
        boundary_integral_ = eval_points_->add_boundary<3>(*q_side_ );
        mesh_ = mesh_full_constructor("{ mesh_file=\"mesh/simplest_cube.msh\", optimize_mesh=false }");
        //std::shared_ptr<DiscreteSpace> ds = std::make_shared<EqualOrderDiscreteSpace>(mesh_, fe_rt0);
        dh_ = std::make_shared<DOFHandlerMultiDim>(*mesh_);
        //dh_->distribute_dofs(ds);
        this->init(eval_points_);

    }

    ~RangesSpeedTest() {
        Profiler::uninitialize();
    }

	void profiler_output() {
		static ofstream os( FilePath("speed_point_ranges_test.log", FilePath::output_file) );
		Profiler::instance()->output(MPI_COMM_WORLD, os);
		os << "" << std::setfill('=') << setw(80) << "" << std::setfill(' ') << endl << endl;
	}

    void add_volume_integral(const DHCellAccessor &cell) {
        unsigned int reg_idx = cell.elm().region_idx().idx();
        unsigned int elm_idx = cell.elm_idx();
        for (uint i=uint(eval_points_->subset_begin(cell.dim(), bulk_integral_->get_subset_idx()));
                  i<uint(eval_points_->subset_end(cell.dim(), bulk_integral_->get_subset_idx())); ++i) {
            this->eval_point_data_.emplace_back(reg_idx, elm_idx, i, cell.local_idx());
        }
        this->eval_point_data_.make_permanent();
    }

    void add_edge_integral(const DHCellSide &cell_side) {
        for( DHCellSide edge_side : cell_side.edge_sides() ) {
            unsigned int reg_idx = edge_side.element().region_idx().idx();
            for (auto p : edge_integral_->points(edge_side, this) ) {
                this->eval_point_data_.emplace_back(reg_idx, edge_side.elem_idx(), p.eval_point_idx(), edge_side.cell().local_idx());
            }
        }
        this->eval_point_data_.make_permanent();
    }

    void add_coupling_integral(const DHCellAccessor &cell, const DHCellSide &ngh_side) {
        unsigned int reg_idx_low = cell.elm().region_idx().idx();
        unsigned int reg_idx_high = ngh_side.element().region_idx().idx();
        for (auto p : coupling_integral_->points(ngh_side, this) ) {
            this->eval_point_data_.emplace_back(reg_idx_high, ngh_side.elem_idx(), p.eval_point_idx(), ngh_side.cell().local_idx());

            auto p_low = p.lower_dim(cell); // equivalent point on low dim cell
           	this->eval_point_data_.emplace_back(reg_idx_low, cell.elm_idx(), p_low.eval_point_idx(), cell.local_idx());
       	}
    }

    void add_boundary_integral(const DHCellSide &bdr_side) {
        unsigned int reg_idx = bdr_side.element().region_idx().idx();
        for (auto p : boundary_integral_->points(bdr_side, this) ) {
            this->eval_point_data_.emplace_back(reg_idx, bdr_side.elem_idx(), p.eval_point_idx(), bdr_side.cell().local_idx());

        	auto p_bdr = p.point_bdr(bdr_side.cond().element_accessor()); // equivalent point on boundary element
        	unsigned int bdr_reg = bdr_side.cond().element_accessor().region_idx().idx();
        	this->eval_point_data_.emplace_back(bdr_reg, bdr_side.cond().bc_ele_idx(), p_bdr.eval_point_idx(), -1);
        }
    }


    Mesh * mesh_;
    std::shared_ptr<DOFHandlerMultiDim> dh_;
    std::shared_ptr<EvalPoints> eval_points_;
    Quadrature *q_bulk_;
    Quadrature *q_side_;
    std::shared_ptr<BulkIntegral> bulk_integral_;
    std::shared_ptr<EdgeIntegral> edge_integral_;
    std::shared_ptr<CouplingIntegral> coupling_integral_;
    std::shared_ptr<BoundaryIntegral> boundary_integral_;
};

TEST_F(RangesSpeedTest, point_ranges) {
	unsigned int counter;
    DHCellAccessor dh_cell(dh_.get(), 3);
    DHCellSide edge_side, ngh_side, bdr_side;
    this->start_elements_update();
    this->add_volume_integral(dh_cell);
    for( DHCellSide cell_side : dh_cell.side_range() ) {
        if ( (cell_side.side().edge().n_sides() == 1) && (cell_side.side().is_boundary()) ) {
            this->add_boundary_integral(cell_side);
            bdr_side = cell_side;
            continue;
        }
        if ( cell_side.n_edge_sides() >= 2 ) {
            this->add_edge_integral(cell_side);
            edge_side = cell_side;
        }
    }
    DHCellAccessor dh_cell_nbr(dh_.get(), 1);
	for( DHCellSide neighb_side : dh_cell_nbr.neighb_sides() ) { // cell -> elm lower dim, neighb_side -> elm higher dim
        if (dh_cell_nbr.dim() != neighb_side.dim()-1) continue;
        this->add_coupling_integral(dh_cell_nbr, neighb_side);
        ngh_side = neighb_side;
        break;
    }
    this->create_patch();
    this->finish_elements_update();

    START_TIMER("TEST_cpp_loop");
    counter = 0;
    for (uint i=0; i<profiler_loops; ++i)
        for (uint j=0; j<dh_cell.dim()+1; ++j) {
            ++counter;
        }
    END_TIMER("TEST_cpp_loop");
    EXPECT_EQ(counter, (dh_cell.dim()+1)*profiler_loops);

    START_TIMER("TEST_bulk_points");
    unsigned int elm_patch_idx = this->position_in_cache(dh_cell.elm_idx());
    counter = 0;
    for (uint i=0; i<profiler_loops; ++i)
        for (auto p : bulk_integral_->points(elm_patch_idx, this)) {
            ++counter;
        }
    END_TIMER("TEST_bulk_points");
    EXPECT_EQ(counter, (dh_cell.dim()+1)*profiler_loops);

    START_TIMER("TEST_edge_points");
    counter = 0;
    for (uint i=0; i<profiler_loops; ++i)
        for (auto p : edge_integral_->points(edge_side, this) ) {
            ++counter;
        }
    END_TIMER("TEST_edge_points");
    EXPECT_EQ(counter, dh_cell.dim()*profiler_loops);

    START_TIMER("TEST_coupling_points");
    counter = 0;
    for (uint i=0; i<profiler_loops; ++i)
        for (auto p : coupling_integral_->points(ngh_side, this) ) {
            ++counter;
        }
    END_TIMER("TEST_coupling_points");
    EXPECT_EQ(counter, (dh_cell.dim())*profiler_loops);

    START_TIMER("TEST_boundary_points");
    counter = 0;
    for (uint i=0; i<profiler_loops; ++i)
        for (auto p : boundary_integral_->points(bdr_side, this) ) {
            ++counter;
        }
    END_TIMER("TEST_boundary_points");
    EXPECT_EQ(counter, (dh_cell.dim())*profiler_loops);

    this->profiler_output();
}


#endif // FLOW123D_RUN_UNIT_BENCHMARKS
