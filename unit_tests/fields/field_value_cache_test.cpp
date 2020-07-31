/*
 * field_value_cache_test.cpp
 *
 *  Created on: Dec 03, 2019
 *      Author: David Flanderka
 *
 *  Tests FieldValueCache and ElementCacheMap classes
 */

#define TEST_USE_MPI
#define FEAL_OVERRIDE_ASSERTS
#include <flow_gtest_mpi.hh>
#include <mesh_constructor.hh>
#include "arma_expect.hh"

#include "fields/eval_points.hh"
#include "fields/eval_subset.hh"
#include "fields/field_value_cache.hh"
#include "fields/field_values.hh"
#include "quadrature/quadrature.hh"
#include "quadrature/quadrature_lib.hh"
#include "fem/dofhandler.hh"
#include "fem/dh_cell_accessor.hh"
#include "mesh/mesh.h"
#include "system/sys_profiler.hh"


typedef FieldValue_<1,1,double> ScalarValue;

class FieldValueCacheTest : public testing::Test, public ElementCacheMap {

public:
    FieldValueCacheTest() {
        FilePath::set_io_dirs(".",UNIT_TESTS_SRC_DIR,"",".");
        Profiler::instance();
        PetscInitialize(0,PETSC_NULL,PETSC_NULL,PETSC_NULL);

        mesh_ = mesh_full_constructor("{mesh_file=\"mesh/cube_2x1.msh\"}");
        dh_ = std::make_shared<DOFHandlerMultiDim>(*mesh_);

        eval_points = std::make_shared<EvalPoints>();
        Quadrature *q_bulk = new QGauss(3, 2);
        Quadrature *q_side = new QGauss(2, 2);
        bulk_eval = eval_points->add_bulk<3>(*q_bulk );
        edge_eval = eval_points->add_edge<3>(*q_side );
        this->init(eval_points);
    }

    ~FieldValueCacheTest() {}

    void add_bulk_points(DHCellAccessor cell) {
        unsigned int reg_idx = cell.elm().region_idx().idx();
        for (auto p : bulk_eval->points(cell, this) ) {
            EvalPointData epd(reg_idx, cell.elm_idx(), p.eval_point_idx());
            this->eval_point_data_.push_back(epd);
        }
        this->eval_point_data_.make_permanent();
        this->add(cell); // temporary: old code, remove
    }

    void add_side_points(DHCellSide cell_side) {
        unsigned int reg_idx = cell_side.element().region_idx().idx();
        for (auto p : edge_eval->points(cell_side, this) ) {
            EvalPointData epd(reg_idx, cell_side.elem_idx(), p.eval_point_idx());
            this->eval_point_data_.push_back(epd);

        	auto p_ghost = p.point_on(cell_side); // point on neghbouring side on one edge
        	unsigned int ghost_reg = p_ghost.dh_cell_side().element().region_idx().idx();
        	EvalPointData epd_ghost(ghost_reg, p_ghost.dh_cell_side().elem_idx(), p_ghost.eval_point_idx());
        	this->eval_point_data_.push_back(epd_ghost);
        }
        this->add(cell_side); // temporary: old code, remove
    }

    Mesh * mesh_;
    std::shared_ptr<DOFHandlerMultiDim> dh_;

    std::shared_ptr<EvalPoints> eval_points;
    std::shared_ptr<BulkIntegral> bulk_eval;
    std::shared_ptr<EdgeIntegral> edge_eval;
};


TEST_F(FieldValueCacheTest, field_value_cache) {
    FieldValueCache<double> value_cache(1, 1);
    unsigned int cache_size = ElementCacheMap::n_cached_elements * eval_points->max_size();
    value_cache.reinit(cache_size);
    value_cache.resize(cache_size);
    EXPECT_EQ(value_cache.size(), eval_points->max_size()*ElementCacheMap::n_cached_elements);

    this->start_elements_update();
    DHCellAccessor dh_cell(dh_.get(), 2);
    this->add(dh_cell);
    for( DHCellSide cell_side : dh_cell.side_range() )
      if ( cell_side.n_edge_sides() >= 2 )
        for( DHCellSide edge_side : cell_side.edge_sides() ) {
            this->add(edge_side);
        }

    this->prepare_elements_to_update();

    // mark used points
    this->mark_used_eval_points(dh_cell, bulk_eval->get_subset_idx(), 4);
    for( DHCellSide cell_side : dh_cell.side_range() )
      if ( cell_side.n_edge_sides() >= 2 )
    	for( DHCellSide edge_side : cell_side.edge_sides() ) {
            this->mark_used_eval_points(edge_side.cell(), edge_eval->get_subset_idx(), 3, 3*edge_side.side_idx());
        }
    this->create_elements_points_map();

    // set value
    unsigned int points_in_cache = update_data_.region_value_cache_range_[update_data_.region_cache_indices_map_.size()];
    EXPECT_EQ(points_in_cache, 16);
    Armor::ArmaMat<double, 1, 1> const_val{0.5};
    for (unsigned int i=0; i<points_in_cache; ++i) value_cache.set(i) = const_val;
    this->finish_elements_update();

    // check value
    dh_cell = this->cache_map_index(dh_cell);
    for(BulkPoint q_point: bulk_eval->points(dh_cell, this)) {
        auto point_val = this->get_value<ScalarValue>(value_cache, dh_cell, q_point.eval_point_idx());
    	EXPECT_DOUBLE_EQ( point_val, const_val(0) );
    }
    for ( DHCellSide cell_side : dh_cell.side_range() )
      if ( cell_side.n_edge_sides() >= 2 )
        for( DHCellSide edge_side : cell_side.edge_sides() )
            for ( EdgePoint q_point : edge_eval->points(edge_side, this) ) {
                auto edge_cell = this->cache_map_index(edge_side.cell());
                auto point_val = this->get_value<ScalarValue>(value_cache, edge_cell, q_point.eval_point_idx());
                EXPECT_DOUBLE_EQ( point_val, const_val(0) );
            }
}


TEST_F(FieldValueCacheTest, element_cache_map) {
    const ElementCacheMap::UpdateCacheHelper &update_cache_data = this->update_cache_data();

    // Test of 2 elements on same region
    this->start_elements_update();
    DHCellAccessor dh_cell1(dh_.get(), 1);
    DHCellAccessor dh_cell2(dh_.get(), 2);
    this->add_bulk_points(dh_cell1);
    this->add_bulk_points(dh_cell2);

    this->prepare_elements_to_update();
    EXPECT_EQ(this->n_elements(), 2);
    EXPECT_EQ(update_cache_data.region_cache_indices_map_.size(), 1);
    EXPECT_TRUE(update_cache_data.region_cache_indices_map_.find(1)!=update_cache_data.region_cache_indices_map_.end());
    EXPECT_EQ(update_cache_data.region_cache_indices_map_.find(1)->second.n_elements_, 2);

    this->create_elements_points_map();
    EXPECT_EQ(update_cache_data.region_cache_indices_map_.size(), 1);
    this->finish_elements_update();
    this->eval_point_data_.reset();

    dh_cell1 = this->cache_map_index(dh_cell1);
    EXPECT_EQ(dh_cell1.element_cache_index(), 0);

    // Test of edge connectivity
    this->start_elements_update();
    //EXPECT_EQ(this->n_elements(), 0);
    EXPECT_EQ(update_cache_data.region_cache_indices_map_.size(), 0);
    for( DHCellSide cell_side : dh_cell2.side_range() )
        if ( cell_side.n_edge_sides() >= 2 )
            for( DHCellSide edge_side : cell_side.edge_sides() ) {
            	this->add_side_points(edge_side);
            }
    this->prepare_elements_to_update();
    //EXPECT_EQ(this->n_elements(), 3); //TODO fix test here
    EXPECT_EQ(update_cache_data.region_cache_indices_map_.size(), 1);
    EXPECT_TRUE(update_cache_data.region_cache_indices_map_.find(1)!=update_cache_data.region_cache_indices_map_.end());
    EXPECT_EQ(update_cache_data.region_cache_indices_map_.find(1)->second.n_elements_, 3);

    for( DHCellSide cell_side : dh_cell2.side_range() )
        if ( cell_side.n_edge_sides() >= 2 )
            for( DHCellSide edge_side : cell_side.edge_sides() ) {
                this->mark_used_eval_points(edge_side.cell(), edge_eval->get_subset_idx(), 3, 3*edge_side.side_idx());
            }
    this->create_elements_points_map();
    EXPECT_EQ(update_cache_data.region_cache_indices_map_.size(), 1);
    EXPECT_EQ(update_cache_data.region_value_cache_range_[0], 0);
    EXPECT_EQ(update_cache_data.region_value_cache_range_[1], 12);
    this->finish_elements_update();
    this->eval_point_data_.reset();
    dh_cell2 = this->cache_map_index(dh_cell2);
    EXPECT_EQ(dh_cell2.element_cache_index(), 1);

    // Test of 3 elements on 2 different regions
    this->start_elements_update();
    DHCellAccessor dh_cell3(dh_.get(), 3);
    DHCellAccessor dh_cell6(dh_.get(), 6);
    this->add_bulk_points(dh_cell1);
    this->add_bulk_points(dh_cell3);
    this->add_bulk_points(dh_cell6);

    this->prepare_elements_to_update();
    EXPECT_EQ(this->n_elements(), 3);
    EXPECT_EQ(update_cache_data.region_cache_indices_map_.size(), 2);
    EXPECT_TRUE(update_cache_data.region_cache_indices_map_.find(1)!=update_cache_data.region_cache_indices_map_.end());
    EXPECT_EQ(update_cache_data.region_cache_indices_map_.find(1)->second.n_elements_, 2);
    EXPECT_EQ(update_cache_data.region_cache_indices_map_.find(3)->second.n_elements_, 1);

    this->create_elements_points_map();
    EXPECT_EQ(update_cache_data.region_cache_indices_map_.size(), 2);
    this->finish_elements_update();
    this->eval_point_data_.reset();
    dh_cell1 = this->cache_map_index(dh_cell1);
    EXPECT_EQ(dh_cell1.element_cache_index(), 1);
}
