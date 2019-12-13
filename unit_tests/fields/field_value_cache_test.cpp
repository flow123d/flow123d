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
#include "fields/field_value_cache.impl.hh"
#include "fields/field_values.hh"
#include "quadrature/quadrature.hh"
#include "quadrature/quadrature_lib.hh"
#include "fem/dofhandler.hh"
#include "fem/dh_cell_accessor.hh"
#include "fem/mapping_p1.hh"
#include "mesh/mesh.h"
#include "system/sys_profiler.hh"
#include "arma_expect.hh"


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


TEST_F(FieldValueCacheTest, field_value_cache) {
    std::shared_ptr<EvalPoints> eval_points = std::make_shared<EvalPoints>();
    Quadrature *q_bulk = new QGauss(3, 2);
    Quadrature *q_side = new QGauss(2, 2);
    std::shared_ptr<EvalSubset> bulk_eval = eval_points->add_bulk<3>(*q_bulk );
    std::shared_ptr<EvalSubset> side_eval = eval_points->add_side<3>(*q_side );

    FieldValueCache<double, double> value_cache(1, 1);
    value_cache.init(eval_points, ElementCacheMap::n_cached_elements);
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
}


TEST_F(FieldValueCacheTest, element_cache_map) {
    ElementCacheMap cache_map;
    cache_map.init(3);
    const ElementCacheMap::UpdateCacheHelper &update_cache_data = cache_map.update_cache_data();

    DHCellAccessor dh_cell1(dh_.get(), 1);
    DHCellAccessor dh_cell2(dh_.get(), 2);
    cache_map.add(dh_cell1);
    cache_map.add(dh_cell2);
    EXPECT_EQ(update_cache_data.added_elements_.size(), 2);

    cache_map.prepare_elements_to_update(mesh_);
    EXPECT_EQ(update_cache_data.preserved_elements_.size(), 0);
    EXPECT_EQ(update_cache_data.region_element_map_.size(), 1);
    EXPECT_TRUE(update_cache_data.region_element_map_.find(1)!=update_cache_data.region_element_map_.end());
    EXPECT_EQ(update_cache_data.region_element_map_.find(1)->second.size(), 2);
    EXPECT_EQ(update_cache_data.region_cache_begin_.find(1)->second, 0);

    cache_map.clear_elements_to_update();
    EXPECT_EQ(update_cache_data.added_elements_.size(), 0);
    EXPECT_EQ(update_cache_data.region_element_map_.size(), 0);

    dh_cell1 = cache_map(dh_cell1);
    EXPECT_EQ(dh_cell1.element_cache_index(), 1);

    DHCellAccessor dh_cell3(dh_.get(), 3);
    DHCellAccessor dh_cell6(dh_.get(), 6);
    cache_map.add(dh_cell1);
    cache_map.add(dh_cell3);
    cache_map.add(dh_cell6);
    EXPECT_EQ(update_cache_data.added_elements_.size(), 3);

    cache_map.prepare_elements_to_update(mesh_);
    EXPECT_EQ(update_cache_data.added_elements_.size(), 2);
    EXPECT_EQ(update_cache_data.preserved_elements_.size(), 1);
    EXPECT_EQ(update_cache_data.region_element_map_.size(), 2);
    EXPECT_TRUE(update_cache_data.region_element_map_.find(1)!=update_cache_data.region_element_map_.end());
    EXPECT_EQ(update_cache_data.region_element_map_.find(1)->second.size(), 1);
    EXPECT_EQ(update_cache_data.region_element_map_.find(3)->second.size(), 1);

    cache_map.clear_elements_to_update();
    dh_cell1 = cache_map(dh_cell1);
    EXPECT_EQ(dh_cell1.element_cache_index(), 0);
}

TEST_F(FieldValueCacheTest, global_coords) {
    std::shared_ptr<EvalPoints> eval_points = std::make_shared<EvalPoints>();
    Quadrature *q_bulk = new QGauss(3, 1);
    Quadrature *q_side = new QGauss(2, 1);
    std::shared_ptr<EvalSubset> bulk_eval = eval_points->add_bulk<3>(*q_bulk );
    std::shared_ptr<EvalSubset> side_eval = eval_points->add_side<3>(*q_side );

    ElementCacheMap cache_map;
    cache_map.init(3);
    const ElementCacheMap::UpdateCacheHelper &update_cache_data = cache_map.update_cache_data();

    { // add cells to cache
        DHCellAccessor dh_cell1(dh_.get(), 1);
        DHCellAccessor dh_cell2(dh_.get(), 2);
        DHCellAccessor dh_cell3(dh_.get(), 3);
        DHCellAccessor dh_cell6(dh_.get(), 6);
        cache_map.add(dh_cell1);
        cache_map.add(dh_cell2);
        cache_map.add(dh_cell3);
        cache_map.add(dh_cell6);
    }

    cache_map.prepare_elements_to_update(mesh_);
    cache_map.clear_elements_to_update();

    MappingP1<3,3> mapping;
    cache_map.compute_global_coords<3>(eval_points, mapping, mesh_);

    DHCellAccessor dh_cell(dh_.get(), 2);
    BulkPoint bulk_point( dh_cell, bulk_eval, eval_points->subset_begin(bulk_eval->get_subset_idx()) );
    EXPECT_ARMA_EQ( cache_map.global_coords(dh_cell.elm_idx(), bulk_point.eval_point_idx()), dh_cell.elm().centre() );
    for (DHCellSide side : dh_cell.side_range()) {
    	auto range_it = side_eval->points(side).begin();
        EXPECT_ARMA_EQ( cache_map.global_coords(side.cell().elm_idx(), range_it->eval_point_idx()), side.side().centre() );
    }

}
