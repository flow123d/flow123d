/*
 * field_const_test.cpp
 *
 *  Created on: Dec 03, 2019
 *      Author: David Flanderka
 *
 *  Tests evaluation of FieldConstant
 */

#define TEST_USE_MPI
#define FEAL_OVERRIDE_ASSERTS
#include <flow_gtest_mpi.hh>
#include <mesh_constructor.hh>
#include "arma_expect.hh"

#include "field_eval_base.hh"
#include "fem/eval_points.hh"
#include "fem/integral_acc.hh"
#include "fem/element_cache_map.hh"
#include "fields/field_values.hh"
#include "fields/field_set.hh"
#include "tools/unit_si.hh"
#include "quadrature/quadrature.hh"
#include "quadrature/quadrature_lib.hh"
#include "fem/dofhandler.hh"
#include "fem/dh_cell_accessor.hh"
#include "mesh/mesh.h"
#include "mesh/accessors.hh"
#include "input/input_type.hh"
#include "input/accessors.hh"
#include "input/reader_to_storage.hh"
#include "system/sys_profiler.hh"

class Balance;


class FieldEvalConstantTest : public FieldEvalBaseTest {
public:

    class EqData : public FieldEvalBaseTest::EqData {
        EqData() : FieldEvalBaseTest::EqData() {}
    };

    FieldEvalConstantTest() : FieldEvalBaseTest() {}

    ~FieldEvalConstantTest() {}
};


TEST_F(FieldEvalConstantTest, full_declaration) {
    string eq_data_input = R"YAML(
    data:
      - region: 3D left
        time: 0.0
        scalar_field: !FieldConstant
          value: 0.5
        vector_field: [1, 2, 3]
        tensor_field: [ [0,1,2], [3,4,5], [6,7,8] ]
      - region: 3D right
        time: 0.0
        scalar_field: !FieldConstant
          value: 15
          unit: "dm"
        vector_field: [4, 5, 6]
        tensor_field: [ [1,2,3], [4,5,6], [7,8,9] ]
    )YAML";

    std::vector<double>       expected_scalar = {0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.5, 1.5, 1.5, 1.5, 1.5, 1.5};
    std::vector<arma::vec3>   expected_vector = {{1, 2, 3}, {1, 2, 3}, {1, 2, 3}, {1, 2, 3}, {1, 2, 3}, {1, 2, 3},
                                                 {4, 5, 6}, {4, 5, 6}, {4, 5, 6}, {4, 5, 6}, {4, 5, 6}, {4, 5, 6}};
    std::vector<arma::mat33>  expected_tensor = {{0.0, 3, 6, 1, 4, 7, 2, 5, 8}, {0.0, 3, 6, 1, 4, 7, 2, 5, 8}, {0.0, 3, 6, 1, 4, 7, 2, 5, 8},
                                                 {0.0, 3, 6, 1, 4, 7, 2, 5, 8}, {0.0, 3, 6, 1, 4, 7, 2, 5, 8}, {0.0, 3, 6, 1, 4, 7, 2, 5, 8},
                                                 {1.0, 4, 7, 2, 5, 8, 3, 6, 9}, {1.0, 4, 7, 2, 5, 8, 3, 6, 9}, {1.0, 4, 7, 2, 5, 8, 3, 6, 9},
                                                 {1.0, 4, 7, 2, 5, 8, 3, 6, 9}, {1.0, 4, 7, 2, 5, 8, 3, 6, 9}, {1.0, 4, 7, 2, 5, 8, 3, 6, 9} };

    this->create_mesh("mesh/cube_2x1.msh");
    this->read_input(eq_data_input);

  	eq_data_->reallocate_cache();

    // BULK fields
  	VecRef<double> ref_scalar(expected_scalar);
   	VecRef<arma::vec3> ref_vector(expected_vector);
   	VecRef<arma::mat33> ref_tensor(expected_tensor);
    EXPECT_TRUE( eval_bulk_field(eq_data_->scalar_field, ref_scalar) );
    EXPECT_TRUE( eval_bulk_field(eq_data_->vector_field, ref_vector) );
    EXPECT_TRUE( eval_bulk_field(eq_data_->tensor_field, ref_tensor) );
}


TEST_F(FieldEvalConstantTest, tensor_single_value) {
    string eq_data_input = R"YAML(
    data:
      - region: 3D left
        time: 0.0
        tensor_field: 3.14
      - region: 3D right
        time: 0.0
        tensor_field: 2.72
    )YAML";

    std::vector<arma::mat33> expected_tensor = {{3.14, 0, 0, 0, 3.14, 0, 0, 0, 3.14}, {3.14, 0, 0, 0, 3.14, 0, 0, 0, 3.14}, {3.14, 0, 0, 0, 3.14, 0, 0, 0, 3.14},
                                                {3.14, 0, 0, 0, 3.14, 0, 0, 0, 3.14}, {3.14, 0, 0, 0, 3.14, 0, 0, 0, 3.14}, {3.14, 0, 0, 0, 3.14, 0, 0, 0, 3.14},
                                                {2.72, 0, 0, 0, 2.72, 0, 0, 0, 2.72}, {2.72, 0, 0, 0, 2.72, 0, 0, 0, 2.72}, {2.72, 0, 0, 0, 2.72, 0, 0, 0, 2.72},
                                                {2.72, 0, 0, 0, 2.72, 0, 0, 0, 2.72}, {2.72, 0, 0, 0, 2.72, 0, 0, 0, 2.72}, {2.72, 0, 0, 0, 2.72, 0, 0, 0, 2.72} };

    this->create_mesh("mesh/cube_2x1.msh");
    this->read_input(eq_data_input);

  	eq_data_->reallocate_cache();

    // BULK fields
   	VecRef<arma::mat33> ref_tensor(expected_tensor);
    EXPECT_TRUE( eval_bulk_field(eq_data_->tensor_field, ref_tensor) );
}


TEST_F(FieldEvalConstantTest, tensor_diagonal_matrix) {
    string eq_data_input = R"YAML(
    data:
      - region: 3D left
        time: 0.0
        tensor_field: [1, 2, 3]
      - region: 3D right
        time: 0.0
        tensor_field: [1.2, 2.3, 3.4]
    )YAML";

    std::vector<arma::mat33> expected_tensor = {{1.0, 0, 0, 0, 2.0, 0, 0, 0, 3.0}, {1.0, 0, 0, 0, 2.0, 0, 0, 0, 3.0}, {1.0, 0, 0, 0, 2.0, 0, 0, 0, 3.0},
                                                {1.0, 0, 0, 0, 2.0, 0, 0, 0, 3.0}, {1.0, 0, 0, 0, 2.0, 0, 0, 0, 3.0}, {1.0, 0, 0, 0, 2.0, 0, 0, 0, 3.0},
                                                {1.2, 0, 0, 0, 2.3, 0, 0, 0, 3.4}, {1.2, 0, 0, 0, 2.3, 0, 0, 0, 3.4}, {1.2, 0, 0, 0, 2.3, 0, 0, 0, 3.4},
                                                {1.2, 0, 0, 0, 2.3, 0, 0, 0, 3.4}, {1.2, 0, 0, 0, 2.3, 0, 0, 0, 3.4}, {1.2, 0, 0, 0, 2.3, 0, 0, 0, 3.4} };

    this->create_mesh("mesh/cube_2x1.msh");
    this->read_input(eq_data_input);

  	eq_data_->reallocate_cache();

    // BULK fields
   	VecRef<arma::mat33> ref_tensor(expected_tensor);
    EXPECT_TRUE( eval_bulk_field(eq_data_->tensor_field, ref_tensor) );
}


TEST_F(FieldEvalConstantTest, tensor_symetric_matrix) {
    string eq_data_input = R"YAML(
    data:
      - region: 3D left
        time: 0.0
        tensor_field: [0.1, 0.2, 0.3, 0.4, 0.5, 0.6]
      - region: 3D right
        time: 0.0
        tensor_field: [2.1, 2.2, 2.3, 2.4, 2.5, 2.6]
    )YAML";

    std::vector<arma::mat33> expected_tensor = {{0.1, 0.2, 0.3, 0.2, 0.4, 0.5, 0.3, 0.5, 0.6}, {0.1, 0.2, 0.3, 0.2, 0.4, 0.5, 0.3, 0.5, 0.6}, {0.1, 0.2, 0.3, 0.2, 0.4, 0.5, 0.3, 0.5, 0.6},
                                                {0.1, 0.2, 0.3, 0.2, 0.4, 0.5, 0.3, 0.5, 0.6}, {0.1, 0.2, 0.3, 0.2, 0.4, 0.5, 0.3, 0.5, 0.6}, {0.1, 0.2, 0.3, 0.2, 0.4, 0.5, 0.3, 0.5, 0.6},
                                                {2.1, 2.2, 2.3, 2.2, 2.4, 2.5, 2.3, 2.5, 2.6}, {2.1, 2.2, 2.3, 2.2, 2.4, 2.5, 2.3, 2.5, 2.6}, {2.1, 2.2, 2.3, 2.2, 2.4, 2.5, 2.3, 2.5, 2.6},
                                                {2.1, 2.2, 2.3, 2.2, 2.4, 2.5, 2.3, 2.5, 2.6}, {2.1, 2.2, 2.3, 2.2, 2.4, 2.5, 2.3, 2.5, 2.6}, {2.1, 2.2, 2.3, 2.2, 2.4, 2.5, 2.3, 2.5, 2.6} };

    this->create_mesh("mesh/cube_2x1.msh");
    this->read_input(eq_data_input);

  	eq_data_->reallocate_cache();

    // BULK fields
   	VecRef<arma::mat33> ref_tensor(expected_tensor);
    EXPECT_TRUE( eval_bulk_field(eq_data_->tensor_field, ref_tensor) );
}

TEST_F(FieldEvalConstantTest, set_time) {
    string eq_data_input = R"YAML(
    data:
      - region: ALL
        time: 0.0
        scalar_field: 0
        bc_scalar_field: 0
      - region: BULK
        time: 2.0
        scalar_field: 1
        bc_scalar_field: 1
      - region: .BOUNDARY
        time: 4.0
        scalar_field: 1
        bc_scalar_field: 0
    )YAML";

    unsigned int        n_times             = 5;   // time steps: { 0, 1, 2, 3, 4 }
    std::vector<double> expected_scalar     = {0.0, 0.0, 1.0, 1.0, 1.0};
    std::vector<bool>   scalar_jump_time    = {true, false, true, false, true};
    std::vector<bool>   bc_scalar_jump_time = {true, false, true, false, true};

    this->create_mesh("mesh/simplest_cube.msh");
    this->read_input(eq_data_input);

    for (unsigned int j=0; j<n_times; j++) {  // time loop
        eq_data_->reallocate_cache();

        // BULK field
        SingleValRef<double> ref_scalar(expected_scalar[j]);
        EXPECT_TRUE( eval_bulk_field(eq_data_->scalar_field, ref_scalar) );
        EXPECT_EQ( scalar_jump_time[j] , eq_data_->scalar_field.is_jump_time() );

        // BOUNDARY field
    	SingleValRef<double> ref_bc_scalar(0.0);
        EXPECT_TRUE( eval_boundary_field(eq_data_->bc_scalar_field, ref_bc_scalar, 3, 0) );
        EXPECT_EQ( bc_scalar_jump_time[j] , eq_data_->bc_scalar_field.is_jump_time() );

        eq_data_->tg_.next_time();
    }
}

TEST_F(FieldEvalConstantTest, copy_constructor) {
    string eq_data_input = R"YAML(
    data:
      - region: BULK
        time: 0.0
        scalar_field: 0
        scalar_ref: 1
      - region: BULK
        time: 1.0
        scalar_ref: 1
      - region: BULK
        time: 2.0
        scalar_field: 1
      - region: BULK
        time: 3.0
        scalar_field: 0
        scalar_ref: 0
    )YAML";

    this->create_mesh("mesh/simplest_cube.msh");
    this->read_input(eq_data_input);

    ScalarField f2(eq_data_->scalar_field);	// default constructor
    eq_data_->scalar_ref = eq_data_->scalar_field; // assignment, should overwrite name "scalar_ref" by name "scalar_field"

    // tg = 0.0
    {
        eq_data_->reallocate_cache();
        f2.set_time(eq_data_->tg_.step(), LimitSide::right);
        SingleValRef<double> ref_scalar(0.0);
        EXPECT_TRUE( eval_bulk_field(f2, ref_scalar) );
        EXPECT_TRUE( eval_bulk_field(eq_data_->scalar_field, ref_scalar) );
    }

    // tg = 1.0
    {
    	eq_data_->tg_.next_time();
        eq_data_->reallocate_cache();
        SingleValRef<double> ref_scalar(0.0);
        EXPECT_TRUE( eval_bulk_field(eq_data_->scalar_field, ref_scalar) );
    }

    // tg = 2.0
	{
        eq_data_->tg_.next_time();
        eq_data_->reallocate_cache();
        SingleValRef<double> ref_scalar(1.0);
        EXPECT_TRUE( eval_bulk_field(eq_data_->scalar_field, ref_scalar) );
        SingleValRef<double> ref_f2(0.0);
        EXPECT_TRUE( eval_bulk_field(f2, ref_f2) );
        EXPECT_TRUE( eval_bulk_field(eq_data_->scalar_ref, ref_scalar) );
	}

}

