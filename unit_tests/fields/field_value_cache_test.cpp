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


class FieldValueCacheTest : public testing::Test {

public:
    FieldValueCacheTest() {
        FilePath::set_io_dirs(".",UNIT_TESTS_SRC_DIR,"",".");
        Profiler::initialize();
        PetscInitialize(0,PETSC_NULL,PETSC_NULL,PETSC_NULL);

        mesh_ = mesh_full_constructor("{mesh_file=\"mesh/cube_2x1.msh\"}");
        dh_ = std::make_shared<DOFHandlerMultiDim>(*mesh_);
    }

    ~FieldValueCacheTest() {}

    Mesh * mesh_;
    std::shared_ptr<DOFHandlerMultiDim> dh_;
};


/*TEST_F(FieldValueCacheTest, field_value_cache) {
    std::shared_ptr<EvalPoints> eval_points = std::make_shared<EvalPoints>();
    Quadrature *q_bulk = new QGauss(3, 2);
    Quadrature *q_side = new QGauss(2, 2);
    std::shared_ptr<EvalSubset> bulk_eval = eval_points->add_bulk_old<3>(*q_bulk );
    std::shared_ptr<EvalSubset> side_eval = eval_points->add_side_old<3>(*q_side );

    FieldValueCache<double, double> value_cache(1, 1);
    value_cache.init(eval_points, bulk_eval->dim(), ElementCacheMap::n_cached_elements);
    EXPECT_EQ(value_cache.n_subsets(), 2);
    EXPECT_EQ(value_cache.subset_begin(0), 0);
    EXPECT_EQ(value_cache.subset_end(0), 4*ElementCacheMap::n_cached_elements);
    EXPECT_EQ(value_cache.subset_size(0), 4*ElementCacheMap::n_cached_elements);
    EXPECT_EQ(value_cache.subset_begin(1), 4*ElementCacheMap::n_cached_elements);
    EXPECT_EQ(value_cache.subset_end(1), 16*ElementCacheMap::n_cached_elements);
    EXPECT_EQ(value_cache.subset_size(1), 12*ElementCacheMap::n_cached_elements);

    value_cache.mark_used(side_eval);
    EXPECT_FALSE(value_cache.used_subsets()[0]);
    EXPECT_TRUE(value_cache.used_subsets()[1]);
}*/


TEST_F(FieldValueCacheTest, element_cache_map) {
    std::shared_ptr<EvalPoints> eval_points = std::make_shared<EvalPoints>();
    Quadrature *q_bulk = new QGauss(3, 2);
    Quadrature *q_side = new QGauss(2, 2);
    std::shared_ptr<BulkIntegral> bulk_eval = eval_points->add_bulk<3>(*q_bulk );
    std::shared_ptr<EdgeIntegral> side_eval = eval_points->add_edge<3>(*q_side );

    ElementCacheMap cache_map;
    cache_map.init(eval_points);
    const ElementCacheMap::UpdateCacheHelper &update_cache_data = cache_map.update_cache_data();

    // Test of 2 elements on same region
    DHCellAccessor dh_cell1(dh_.get(), 1);
    DHCellAccessor dh_cell2(dh_.get(), 2);
    cache_map.add(dh_cell1);
    cache_map.add(dh_cell2);
    EXPECT_EQ(update_cache_data.added_elements_.size(), 2);

    cache_map.prepare_elements_to_update(mesh_);
    EXPECT_EQ(update_cache_data.region_element_map_.size(), 1);
    EXPECT_TRUE(update_cache_data.region_element_map_.find(1)!=update_cache_data.region_element_map_.end());
    EXPECT_EQ(update_cache_data.region_element_map_.find(1)->second.size(), 2);
    EXPECT_EQ(update_cache_data.region_cache_begin_.find(1)->second, 0);

    cache_map.clear_elements_to_update();
    EXPECT_EQ(update_cache_data.added_elements_.size(), 0);
    EXPECT_EQ(update_cache_data.region_element_map_.size(), 0);

    dh_cell1 = cache_map(dh_cell1);
    EXPECT_EQ(dh_cell1.element_cache_index(), 1);

    // Test of edge connectivity
    for( DHCellSide cell_side : dh_cell2.side_range() )
        if ( cell_side.n_edge_sides() >= 2 )
            for( DHCellSide edge_side : cell_side.edge_sides() ) {
                cache_map.add(edge_side);
            }
    EXPECT_EQ(update_cache_data.added_elements_.size(), 3);
    cache_map.prepare_elements_to_update(mesh_);
    EXPECT_EQ(update_cache_data.region_element_map_.size(), 1);
    EXPECT_TRUE(update_cache_data.region_element_map_.find(1)!=update_cache_data.region_element_map_.end());
    EXPECT_EQ(update_cache_data.region_element_map_.find(1)->second.size(), 3);

    for( DHCellSide cell_side : dh_cell2.side_range() )
        if ( cell_side.n_edge_sides() >= 2 )
            for( DHCellSide edge_side : cell_side.edge_sides() ) {
    	        cache_map.mark_used_eval_points(edge_side.cell(), side_eval->get_subset_idx(), 3, 3*edge_side.side_idx());
            }
    cache_map.clear_elements_to_update();
    dh_cell1 = cache_map(dh_cell2);
    EXPECT_EQ(dh_cell2.element_cache_index(), 1);

    // Test of 3 elements on 2 different regions
    DHCellAccessor dh_cell3(dh_.get(), 3);
    DHCellAccessor dh_cell6(dh_.get(), 6);
    cache_map.add(dh_cell1);
    cache_map.add(dh_cell3);
    cache_map.add(dh_cell6);
    EXPECT_EQ(update_cache_data.added_elements_.size(), 3);

    cache_map.prepare_elements_to_update(mesh_);
    EXPECT_EQ(update_cache_data.region_element_map_.size(), 2);
    EXPECT_TRUE(update_cache_data.region_element_map_.find(1)!=update_cache_data.region_element_map_.end());
    EXPECT_EQ(update_cache_data.region_element_map_.find(1)->second.size(), 2);
    EXPECT_EQ(update_cache_data.region_element_map_.find(3)->second.size(), 1);

    cache_map.clear_elements_to_update();
    dh_cell1 = cache_map(dh_cell1);
    EXPECT_EQ(dh_cell1.element_cache_index(), 1);
}
