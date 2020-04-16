/*
 * field_model_test.cpp
 *
 *  Created on: Feb 26, 2020
 *      Author: David Flanderka
 *
 *  Tests FieldModel class and expand tuple to arguments.
 */

#define FEAL_OVERRIDE_ASSERTS

#include <flow_gtest.hh>
#include <mesh_constructor.hh>
#include "arma_expect.hh"
#include <iostream>
#include <tuple>
#include <string>
#include <utility>
#include <type_traits>

#include "fields/field_model.hh"
#include "mesh/accessors.hh"
#include "quadrature/quadrature.hh"
#include "quadrature/quadrature_lib.hh"


typedef Field<3, FieldValue<3>::Scalar > ScalarField;
typedef Field<3, FieldValue<3>::VectorFixed > VectorField;


using Scalar = typename arma::Col<double>::template fixed<1>;
using Vector = arma::vec3;
using Tensor = arma::mat33;


// Functor with resolution 'scalar * vector'
Vector fn_product(Scalar a, Vector v) {
    return a(0) * v;
}


// Functor with resolution 'vector + scalar * vector'
Vector fn_other(Vector a, Scalar c, Vector b) {
    return a + c(0) * b;
}


// Test of FieldModel - used objects and functionalities in field_model.hh.
TEST(FieldModelTest, own_model) {
    unsigned int n_items = 10; // number of tested items

    ScalarField f_scal;
    VectorField f_vec;
    FieldValueCache<typename FieldValue<3>::VectorFixed::element_type> fvc(FieldValue<3>::VectorFixed::NRows_, FieldValue<3>::VectorFixed::NCols_);
    ElementCacheMap elm_cache_map;
    std::vector< ElementAccessor<3> > element_set;

    // initialize field caches
    std::shared_ptr<EvalPoints> eval_points = std::make_shared<EvalPoints>();
    Quadrature *q_bulk = new QGauss(3, 2);
    Quadrature *q_side = new QGauss(2, 2);
    eval_points->add_bulk<3>(*q_bulk );
    eval_points->add_edge<3>(*q_side );
    elm_cache_map.init(eval_points);
    fvc.init(eval_points, ElementCacheMap::n_cached_elements);
    f_scal.cache_allocate(eval_points);
    f_vec.cache_allocate(eval_points);

    // fill field caches
    elm_cache_map.start_elements_update();
    auto &cache_data = elm_cache_map.update_cache_data();
    cache_data.region_cache_indices_map_.insert( {1, ElementCacheMap::RegionData()} );
    cache_data.region_cache_indices_map_.find(1)->second.pos_ = 0;
    cache_data.region_value_cache_range_[0] = 0;
    cache_data.region_value_cache_range_[1] = 10;
    arma::mat::fixed<1,1> scalar_val;
    arma::mat::fixed<3,1> vector_val;
    for (unsigned int i=0; i<n_items; ++i) {
        scalar_val(0,0) = 1.0 + i*0.5;
        vector_val(0,0) = 1.5 + 2*i;
        vector_val(1,0) = i + 0.1;
        vector_val(2,0) = 0.5 + i%2;
        f_scal.value_cache().data().set(i) = scalar_val;
        f_vec.value_cache().data().set(i) = vector_val;
    }

    {
        // FieldModel scalar * vector
        std::vector<arma::vec3> expected_vals = {{  1.50,  0.10, 0.50},
                                                 {  5.25,  1.65, 2.25},
                                                 { 11.00,  4.20, 1.00},
                                                 { 18.75,  7.75, 3.75},
                                                 { 28.50, 12.30, 1.50},
                                                 { 40.25, 17.85, 5.25},
                                                 { 54.00, 24.40, 2.00},
                                                 { 69.75, 31.95, 6.75},
                                                 { 87.50, 40.50, 2.50},
                                                 {107.25, 50.05, 8.25}};

    	auto f_product = Model<3, FieldValue<3>::VectorFixed>::create(fn_product, f_scal, f_vec);
        f_product->cache_update(fvc, elm_cache_map, 1);
        for (unsigned int i=0; i<n_items; ++i) {
            auto val = fvc.data().template mat<3, 1>(i);
            EXPECT_ARMA_EQ(val, expected_vals[i]);
        }
    }

    {
        // FieldModel vector + scalar * vector
        std::vector<arma::vec3> expected_vals = {{  3.00,  0.20, 1.00},
                                                 {  8.75,  2.75, 3.75},
                                                 { 16.50,  6.30, 1.50},
                                                 { 26.25, 10.85, 5.25},
                                                 { 38.00, 16.40, 2.00},
                                                 { 51.75, 22.95, 6.75},
                                                 { 67.50, 30.50, 2.50},
                                                 { 85.25, 39.05, 8.25},
                                                 {105.00, 48.60, 3.00},
                                                 {126.75, 59.15, 9.75}};

        auto f_other = Model<3, FieldValue<3>::VectorFixed>::create(fn_other, f_vec, f_scal, f_vec);
        f_other->cache_update(fvc, elm_cache_map, 1);
        for (unsigned int i=0; i<n_items; ++i) {
            auto val = fvc.data().template mat<3, 1>(i);
            EXPECT_ARMA_EQ(val, expected_vals[i]);
        }
    }

}
