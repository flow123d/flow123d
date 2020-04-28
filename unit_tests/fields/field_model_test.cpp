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
#include "fields/multi_field.hh"
#include "mesh/accessors.hh"
#include "mesh/mesh.h"
#include "quadrature/quadrature.hh"
#include "quadrature/quadrature_lib.hh"


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


// Test of FieldModel - simple test without MultiFields (static method Model::create)
TEST(FieldModelTest, create) {
	FilePath::set_io_dirs(".",UNIT_TESTS_SRC_DIR,"",".");
	PetscInitialize(0,PETSC_NULL,PETSC_NULL,PETSC_NULL);
    unsigned int n_items = 10; // number of tested items

    Field<3, FieldValue<3>::Scalar > f_scal;
    Field<3, FieldValue<3>::VectorFixed > f_vec;
    ElementCacheMap elm_cache_map;

    // initialize field caches
    std::shared_ptr<EvalPoints> eval_points = std::make_shared<EvalPoints>();
    Quadrature *q_bulk = new QGauss(3, 2);
    Quadrature *q_side = new QGauss(2, 2);
    eval_points->add_bulk<3>(*q_bulk );
    eval_points->add_edge<3>(*q_side );
    elm_cache_map.init(eval_points);
    f_scal.cache_allocate(eval_points);
    f_vec.cache_allocate(eval_points);

    // Create FieldModel (descendant of FieladAlgoBase) set to Field
    Mesh *mesh = mesh_full_constructor("{mesh_file=\"mesh/cube_2x1.msh\"}");
    TimeGovernor tg(0.0, 1.0);
    auto f_product_ptr = Model<3, FieldValue<3>::VectorFixed>::create(fn_product, f_scal, f_vec);
    Field<3, FieldValue<3>::VectorFixed > f_product;
    f_product.set_mesh( *mesh );
    f_product.set_field(mesh->region_db().get_region_set("ALL"), f_product_ptr);
    f_product.cache_allocate(eval_points);
    f_product.set_time(tg.step(), LimitSide::right);
    // Same as previous but with other functor
    auto f_other_ptr = Model<3, FieldValue<3>::VectorFixed>::create(fn_other, f_vec, f_scal, f_vec);
    Field<3, FieldValue<3>::VectorFixed > f_other;
    f_other.set_mesh( *mesh );
    f_other.set_field(mesh->region_db().get_region_set("ALL"), f_other_ptr);
    f_other.cache_allocate(eval_points);
    f_other.set_time(tg.step(), LimitSide::right);

    // fill field caches
    elm_cache_map.start_elements_update();
    auto &cache_data = elm_cache_map.update_cache_data();
    cache_data.region_cache_indices_map_.insert( {1, ElementCacheMap::RegionData()} );
    cache_data.region_cache_indices_range_.insert( {1, 0} );
    cache_data.region_cache_indices_range_.find(1)->second = 0;
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

        f_product.cache_update(elm_cache_map);
        for (unsigned int i=0; i<n_items; ++i) {
            auto val = f_product.value_cache().data().template mat<3, 1>(i);
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

        f_other.cache_update(elm_cache_map);
        for (unsigned int i=0; i<n_items; ++i) {
            auto val = f_other.value_cache().data().template mat<3, 1>(i);
            EXPECT_ARMA_EQ(val, expected_vals[i]);
        }
    }

}



typedef MultiField<3, FieldValue<3>::Scalar> MultiField3Comp;
typedef MultiField3Comp::SubFieldBaseType ScalarField;


// Functor with resolution 'scalar * multi'
Scalar multi_product(Scalar a, Scalar v) {
    return a * v;
}


// Functor with resolution 'multi + scalar * multi'
Scalar multi_other(Scalar a, Scalar c, Scalar b) {
    return a + c * b;
}


// Test of FieldModel - test of MultiFields (static method Model::create_multi)
TEST(FieldModelTest, create_multi) {
    typedef FieldAlgorithmBase<3, FieldValue<3>::Scalar> FieldBaseType;
    typedef std::shared_ptr< FieldBaseType > FieldBasePtr;

    FilePath::set_io_dirs(".",UNIT_TESTS_SRC_DIR,"",".");
	PetscInitialize(0,PETSC_NULL,PETSC_NULL,PETSC_NULL);
    unsigned int n_items = 10; // number of tested items

    Field<3, FieldValue<3>::Scalar> f_scal;
    MultiField3Comp f_multi;
    ElementCacheMap elm_cache_map;
    std::vector<string> component_names = { "comp_0", "comp_1", "comp_2" };
    Mesh *mesh = mesh_full_constructor("{mesh_file=\"mesh/cube_2x1.msh\"}");
    TimeGovernor tg(0.0, 1.0);

    // initialize field caches
    std::shared_ptr<EvalPoints> eval_points = std::make_shared<EvalPoints>();
    Quadrature *q_bulk = new QGauss(3, 2);
    Quadrature *q_side = new QGauss(2, 2);
    eval_points->add_bulk<3>(*q_bulk );
    eval_points->add_edge<3>(*q_side );
    elm_cache_map.init(eval_points);
    f_scal.name("field_scalar");
    f_scal.cache_allocate(eval_points);
    f_multi.name("field_multi");
    f_multi.set_components(component_names);
    f_multi.set_mesh( *mesh );
    std::vector<FieldBasePtr> field_vec;
    for (uint i=0; i<3; ++i) {
        field_vec.push_back( std::make_shared< FieldConstant<3, FieldValue<3>::Scalar> >() );
    }
    f_multi.set_fields(mesh->region_db().get_region_set("ALL"), field_vec);
    f_multi.cache_allocate(eval_points); // cache_allocate must be called after set_fields!!

    // Create FieldModel (descendant of FieladAlgoBase) set to Field
    auto f_product_ptr = Model<3, FieldValue<3>::Scalar>::create_multi(multi_product, f_scal, f_multi);
    MultiField<3, FieldValue<3>::Scalar> f_product;
    f_product.set_components(component_names);
    f_product.set_mesh( *mesh );
    f_product.set_fields(mesh->region_db().get_region_set("ALL"), f_product_ptr);
    f_product.cache_allocate(eval_points);
    f_product.set_time(tg.step(), LimitSide::right);
    // Same as previous but with other functor
    auto f_other_ptr = Model<3, FieldValue<3>::Scalar>::create_multi(multi_other, f_multi, f_scal, f_multi);
    MultiField<3, FieldValue<3>::Scalar> f_other;
    f_other.set_components(component_names);
    f_other.set_mesh( *mesh );
    f_other.set_fields(mesh->region_db().get_region_set("ALL"), f_other_ptr);
    f_other.cache_allocate(eval_points);
    f_other.set_time(tg.step(), LimitSide::right);

    // fill field caches
    elm_cache_map.start_elements_update();
    auto &cache_data = elm_cache_map.update_cache_data();
    cache_data.region_cache_indices_map_.insert( {1, ElementCacheMap::RegionData()} );
    cache_data.region_cache_indices_range_.insert( {1, 0} );
    cache_data.region_cache_indices_range_.find(1)->second = 0;
    cache_data.region_value_cache_range_[0] = 0;
    cache_data.region_value_cache_range_[1] = 10;
    arma::mat::fixed<1,1> scalar_val;
    std::vector< arma::mat::fixed<1,1> > multi_val;
    multi_val.resize( f_multi.size() );
    for (unsigned int i=0; i<n_items; ++i) {
        scalar_val(0,0) = 1.0 + i*0.5;
        multi_val[0](0,0) = 1.5 + 2*i;
        multi_val[1](0,0) = i + 0.1;
        multi_val[2](0,0) = 0.5 + i%2;
        f_scal.value_cache().data().set(i) = scalar_val;
        for (unsigned int j=0; j<f_multi.size(); ++j)
            f_multi[j].value_cache().data().set(i) = multi_val[j];
    }

    {
        // FieldModel scalar * multi
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

        f_product.cache_update(elm_cache_map);
        for (unsigned int i_cache=0; i_cache<n_items; ++i_cache) {
            for (unsigned int i_subfield=0; i_subfield<f_product.size(); ++i_subfield) {
                auto val = f_product[i_subfield].value_cache().data().template mat<1, 1>(i_cache);
                EXPECT_DOUBLE_EQ(expected_vals[i_cache](i_subfield), val(0));
            }
        }
    }

    {
        // FieldModel multi + scalar * multi
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

        f_other.cache_update(elm_cache_map);
        for (unsigned int i_cache=0; i_cache<n_items; ++i_cache) {
            for (unsigned int i_subfield=0; i_subfield<f_other.size(); ++i_subfield) {
                auto val = f_other[i_subfield].value_cache().data().template mat<1, 1>(i_cache);
                EXPECT_DOUBLE_EQ(expected_vals[i_cache](i_subfield), val(0));
            }
        }
    }

}
