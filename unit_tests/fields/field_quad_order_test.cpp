/*
 * field_quad_order_test.cpp
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
#include "fem/eval_subset.hh"
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


class FieldEvalTest : public FieldEvalBaseTest {
public:

    FieldEvalTest() : FieldEvalBaseTest(2, 0) {}

    ~FieldEvalTest() {}
};


TEST_F(FieldEvalTest, full_declaration) {
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

