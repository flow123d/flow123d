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
#include "fields/field_constant.hh"
#include "fem/integral_acc.hh"
#include "mesh/accessors.hh"
#include "mesh/mesh.h"
#include "quadrature/quadrature.hh"
#include "quadrature/quadrature_lib.hh"
#include "system/sys_profiler.hh"


using Scalar = double;
using Vector = arma::vec3;
using Tensor = arma::mat33;

class FieldModelTest : public testing::Test, public ElementCacheMap {
public:
    virtual void SetUp() {
        Profiler::instance();
    	FilePath::set_io_dirs(".",UNIT_TESTS_SRC_DIR,"",".");
    	PetscInitialize(0,PETSC_NULL,PETSC_NULL,PETSC_NULL);
        n_items = 10; // number of tested items

        mesh = mesh_full_constructor("{ mesh_file=\"mesh/cube_2x1.msh\", optimize_mesh=false }");

        expected_product = {{  1.50,  0.10, 0.50},
                            {  5.25,  1.65, 2.25},
                            { 11.00,  4.20, 1.00},
                            { 18.75,  7.75, 3.75},
                            { 28.50, 12.30, 1.50},
                            { 40.25, 17.85, 5.25},
                            { 54.00, 24.40, 2.00},
                            { 69.75, 31.95, 6.75},
                            { 87.50, 40.50, 2.50},
                            {107.25, 50.05, 8.25}};
        expected_other =   {{  3.00,  0.20, 1.00},
                            {  8.75,  2.75, 3.75},
                            { 16.50,  6.30, 1.50},
                            { 26.25, 10.85, 5.25},
                            { 38.00, 16.40, 2.00},
                            { 51.75, 22.95, 6.75},
                            { 67.50, 30.50, 2.50},
                            { 85.25, 39.05, 8.25},
                            {105.00, 48.60, 3.00},
                            {126.75, 59.15, 9.75}};
    }

    virtual void TearDown() {
    	if (mesh != nullptr) delete mesh;
        Profiler::uninitialize();
    }

    void init_field_caches() {
        eval_points = std::make_shared<EvalPoints>();
        q_bulk = new QGauss(3, 2);
        q_side = new QGauss(2, 2);
        std::shared_ptr<BulkIntegral> mass_eval = std::make_shared<BulkIntegral>(eval_points, q_bulk, 3);
        std::shared_ptr<EdgeIntegral> side_eval = std::make_shared<EdgeIntegral>(eval_points, q_side, 3);
        this->init(eval_points);
    }

    void fill_cache_data() {
    	this->regions_starts_.emplace_back(0);
    	this->regions_starts_.emplace_back(1);
    	this->regions_starts_.make_permanent();
    	this->element_starts_.emplace_back(0);
    	this->element_starts_.emplace_back(n_items);
    	this->element_starts_.make_permanent();
    	for (unsigned int i=0; i<n_items; ++i) {
    	    this->eval_point_data_.emplace_back(1, 1, i, 0);
    	}
    	this->eval_point_data_.make_permanent();
    }

    unsigned int n_items; // number of tested items
    Mesh *mesh;
    std::shared_ptr<EvalPoints> eval_points;
    Quadrature *q_bulk;
    Quadrature *q_side;
    std::vector<arma::vec3> expected_product; // Expected values of f_product
    std::vector<arma::vec3> expected_other; // Expected values of f_other
};

// Functor with resolution 'scalar * vector'
Vector fn_product(Scalar a, Vector v) {
    return a * v;
}


// Functor with resolution 'vector + scalar * vector'
Vector fn_other(Vector a, Scalar c, Vector b) {
    return a + c * b;
}


// Test of FieldModel - simple test without MultiFields (static method Model::create)
TEST_F(FieldModelTest, create) {
    Field<3, FieldValue<3>::Scalar > f_scal;
    Field<3, FieldValue<3>::VectorFixed > f_vec;
    TimeGovernor tg(0.0, 1.0);

    // initialize field caches
    this->init_field_caches();

    // Create FieldModel (descendant of FieladAlgoBase) set to Field
    auto f_product_ptr = Model<3, FieldValue<3>::VectorFixed>::create(fn_product, f_scal, f_vec);
    Field<3, FieldValue<3>::VectorFixed > f_product;
    f_product.set_mesh( *mesh );
    f_product.set(f_product_ptr, 0.0);
    f_product.set_time(tg.step(), LimitSide::right);
    // Same as previous but with other functor
    auto f_other_ptr = Model<3, FieldValue<3>::VectorFixed>::create(fn_other, f_vec, f_scal, f_vec);
    Field<3, FieldValue<3>::VectorFixed > f_other;
    f_other.set_mesh( *mesh );
    f_other.set(f_other_ptr, 0.0);
    f_other.set_time(tg.step(), LimitSide::right);

    // fill field caches
    this->start_elements_update();
    this->fill_cache_data();
    arma::mat::fixed<3,1> vector_val;
    for (unsigned int i=0; i<n_items; ++i) {
        f_scal.value_cache()->set(i) = 1.0 + i*0.5;
        vector_val(0,0) = 1.5 + 2*i;
        vector_val(1,0) = i + 0.1;
        vector_val(2,0) = 0.5 + i%2;
        f_vec.value_cache()->set(i) = vector_val;
    }

    {
        // FieldModel scalar * vector
        f_product.cache_update(*this, 0);
        for (unsigned int i=0; i<n_items; ++i) {
            auto val = f_product.value_cache()->template mat<3, 1>(i);
            EXPECT_ARMA_EQ(val, expected_product[i]);
        }
    }

    {
        // FieldModel vector + scalar * vector
        f_other.cache_update(*this, 0);
        for (unsigned int i=0; i<n_items; ++i) {
            auto val = f_other.value_cache()->template mat<3, 1>(i);
            EXPECT_ARMA_EQ(val, expected_other[i]);
        }
    }

}




// Functor with resolution 'scalar * multi'
Scalar multi_product(Scalar a, Scalar v) {
    return a * v;
}


// Functor with resolution 'multi + scalar * multi'
Scalar multi_other(Scalar a, Scalar c, Scalar b) {
    return a + c * b;
}

// add test multi<vector>


// Test of FieldModel - test of scalar MultiFields (static method Model::create_multi)
TEST_F(FieldModelTest, create_multi_scalar) {
    typedef FieldAlgorithmBase<3, FieldValue<3>::Scalar> FieldBaseType;
    typedef std::shared_ptr< FieldBaseType > FieldBasePtr;

    Field<3, FieldValue<3>::Scalar> f_scal;
    MultiField<3, FieldValue<3>::Scalar> f_multi;
    std::vector<string> component_names = { "comp_0", "comp_1", "comp_2" };
    TimeGovernor tg(0.0, 1.0);

    // initialize field caches
    this->init_field_caches();
    f_scal.name("field_scalar");
    f_scal.set_components(component_names);
    f_multi.name("field_multi");
    f_multi.set_components(component_names);
    f_multi.set_mesh( *mesh );
    std::vector<FieldBasePtr> field_vec;
    for (uint i=0; i<3; ++i) {
        field_vec.push_back( std::make_shared< FieldConstant<3, FieldValue<3>::Scalar> >() );
    }
    f_multi.set(field_vec, 0.0);

    // Create FieldModel (descendant of FieladAlgoBase) set to Field
    auto f_product_ptr = Model<3, FieldValue<3>::Scalar>::create_multi(multi_product, f_scal, f_multi);
    MultiField<3, FieldValue<3>::Scalar> f_product;
    f_product.set_components(component_names);
    f_product.set_mesh( *mesh );
    f_product.set(f_product_ptr, 0.0);
    f_product.set_time(tg.step(), LimitSide::right);
    // Same as previous but with other functor
    auto f_other_ptr = Model<3, FieldValue<3>::Scalar>::create_multi(multi_other, f_multi, f_scal, f_multi);
    MultiField<3, FieldValue<3>::Scalar> f_other;
    f_other.set_components(component_names);
    f_other.set_mesh( *mesh );
    f_other.set(f_other_ptr, 0.0);
    f_other.set_time(tg.step(), LimitSide::right);

    // fill field caches
    this->start_elements_update();
    this->fill_cache_data();
    std::vector< arma::mat::fixed<1,1> > multi_val;
    multi_val.resize( f_multi.size() );
    for (unsigned int i=0; i<n_items; ++i) {
        f_scal.value_cache()->set(i) = 1.0 + i*0.5;
        multi_val[0](0,0) = 1.5 + 2*i;
        multi_val[1](0,0) = i + 0.1;
        multi_val[2](0,0) = 0.5 + i%2;
        for (unsigned int j=0; j<f_multi.size(); ++j)
            f_multi[j].value_cache()->set(i) = multi_val[j];
    }

    {
        // FieldModel scalar * multi
    	for (unsigned int i=0; i<f_product.size(); ++i)
    	    f_product[i].cache_update(*this, 0);
        for (unsigned int i_cache=0; i_cache<n_items; ++i_cache) {
            for (unsigned int i_subfield=0; i_subfield<f_product.size(); ++i_subfield) {
                auto val = f_product[i_subfield].value_cache()->template mat<1, 1>(i_cache);
                EXPECT_DOUBLE_EQ(expected_product[i_cache](i_subfield), val(0));
            }
        }
    }

    {
        // FieldModel multi + scalar * multi
    	for (unsigned int i=0; i<f_other.size(); ++i)
    	    f_other[i].cache_update(*this, 0);
        for (unsigned int i_cache=0; i_cache<n_items; ++i_cache) {
            for (unsigned int i_subfield=0; i_subfield<f_other.size(); ++i_subfield) {
                auto val = f_other[i_subfield].value_cache()->template mat<1, 1>(i_cache);
                EXPECT_DOUBLE_EQ(expected_other[i_cache](i_subfield), val(0));
            }
        }
    }

}

// Test of FieldModel - test of vector MultiFields (static method Model::create_multi)
TEST_F(FieldModelTest, create_multi_vector) {
    typedef FieldAlgorithmBase<3, FieldValue<3>::VectorFixed> FieldBaseType;
    typedef std::shared_ptr< FieldBaseType > FieldBasePtr;

    Field<3, FieldValue<3>::Scalar> f_scal;
    MultiField<3, FieldValue<3>::VectorFixed> f_multi;
    std::vector<string> component_names = { "comp_0", "comp_1", "comp_2" };
    TimeGovernor tg(0.0, 1.0);

    // initialize field caches
    this->init_field_caches();
    f_scal.name("field_scalar");
    f_scal.set_components(component_names);
    f_multi.name("field_multi");
    f_multi.set_components(component_names);
    f_multi.set_mesh( *mesh );
    std::vector<FieldBasePtr> field_vec;
    for (uint i=0; i<3; ++i) {
        field_vec.push_back( std::make_shared< FieldConstant<3, FieldValue<3>::VectorFixed> >() );
    }
    f_multi.set(field_vec, 0.0);

    // Create FieldModel (descendant of FieladAlgoBase) set to Field
    auto f_product_ptr = Model<3, FieldValue<3>::VectorFixed>::create_multi(fn_product, f_scal, f_multi);
    MultiField<3, FieldValue<3>::VectorFixed> f_product;
    f_product.set_components(component_names);
    f_product.set_mesh( *mesh );
    f_product.set(f_product_ptr, 0.0);
    f_product.set_time(tg.step(), LimitSide::right);
    // Same as previous but with other functor
    auto f_other_ptr = Model<3, FieldValue<3>::VectorFixed>::create_multi(fn_other, f_multi, f_scal, f_multi);
    MultiField<3, FieldValue<3>::VectorFixed> f_other;
    f_other.set_components(component_names);
    f_other.set_mesh( *mesh );
    f_other.set(f_other_ptr, 0.0);
    f_other.set_time(tg.step(), LimitSide::right);

    // fill field caches
    this->start_elements_update();
    this->fill_cache_data();
    arma::mat::fixed<3,1> vector_val;
    for (unsigned int i=0; i<n_items; ++i) {
        f_scal.value_cache()->set(i) = 1.0 + i*0.5;
        vector_val(0,0) = 1.5 + 2*i;
        vector_val(1,0) = i + 0.1;
        vector_val(2,0) = 0.5 + i%2;
        for (unsigned int j=0; j<f_multi.size(); ++j)
            f_multi[j].value_cache()->set(i) = vector_val;
    }

    {
        // FieldModel scalar * multi
    	for (unsigned int i=0; i<f_product.size(); ++i)
            f_product[i].cache_update(*this, 0);
        for (unsigned int i_cache=0; i_cache<n_items; ++i_cache) {
            for (unsigned int i_subfield=0; i_subfield<f_product.size(); ++i_subfield) {
                auto val = f_product[i_subfield].value_cache()->template mat<3, 1>(i_cache);
                EXPECT_ARMA_EQ(val, expected_product[i_cache]);
            }
        }
    }

    {
        // FieldModel multi + scalar * multi
    	for (unsigned int i=0; i<f_other.size(); ++i)
    	    f_other[i].cache_update(*this, 0);
        for (unsigned int i_cache=0; i_cache<n_items; ++i_cache) {
            for (unsigned int i_subfield=0; i_subfield<f_other.size(); ++i_subfield) {
                auto val = f_other[i_subfield].value_cache()->template mat<3, 1>(i_cache);
                EXPECT_ARMA_EQ(val, expected_other[i_cache]);
            }
        }
    }

}
