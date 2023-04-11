/*
 * field_time_function_test.cpp
 *
 *  Created on: Dec 18, 2012
 *      Author: jb
 */



#define TEST_USE_MPI
#define FEAL_OVERRIDE_ASSERTS
#include <flow_gtest_mpi.hh>
#include <mesh_constructor.hh>
#include "arma_expect.hh"


#include "field_eval_base.hh"
#include "fields/eval_points.hh"
#include "fields/eval_subset.hh"
#include "fields/field_value_cache.hh"
#include "fields/field_values.hh"
#include "fields/field_set.hh"
#include "fields/field_constant.hh"
#include "fields/field_time_function.hh"
#include "fields/table_function.hh"
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
#include "coupling/generic_assembly.hh"
#include "coupling/assembly_base.hh"

class Balance;


class FieldTimeFunctionTest : public FieldEvalBaseTest {
public:

    class EqData : public FieldEvalBaseTest::EqData {
        EqData() : FieldEvalBaseTest::EqData() {}
    };

    FieldTimeFunctionTest() : FieldEvalBaseTest() {}

    ~FieldTimeFunctionTest() {}
};


TEST_F(FieldTimeFunctionTest, evaluate) {
    string eq_data_input = R"YAML(
    data:
      - region: BULK
        time: 0.0
        scalar_field: !FieldTimeFunction
          time_function:
            - - 0.0
              - 0.5
            - - 2.0
              - 1.0
            - - 4.0
              - 3.0
        vector_field: !FieldTimeFunction
          time_function:
            - - 0.0
              - [0.5, 1.5, 2.0]
            - - 2.0
              - [1.0, 1.5, 1.0]
            - - 4.0
              - [3.0, 1.5, 5.0]
        tensor_field: !FieldTimeFunction
          time_function:
            - - 0.0
              - [ [1,3,4], [0,3,4], [1,6,6] ]
            - - 2.0
              - [ [3,3,6], [0,2,5], [2,7,4] ]
            - - 4.0
              - [ [5,3,4], [2,3,7], [5,6,3] ]
    )YAML";

    unsigned int             n_times         = 6;   // time steps: { 0, 1, 2, 3, 4, 5 }
    std::vector<double>      expected_scalar = {0.5, 0.75, 1.0, 2.0, 3.0, 3.0};
    std::vector<arma::vec3>  expected_vector = {{0.5, 1.5, 2.0}, {0.75, 1.5, 1.5}, {1.0, 1.5, 1.0}, {2.0, 1.5, 3.0}, {3.0, 1.5, 5.0}, {3.0, 1.5, 5.0}};
    std::vector<arma::mat33> expected_tensor = {{1.0, 0.0, 1.0, 3.0, 3.0, 6.0, 4.0, 4.0, 6.0}, {2.0, 0.0, 1.5, 3.0, 2.5, 6.5, 5.0, 4.5, 5.0},
                                                {3.0, 0.0, 2.0, 3.0, 2.0, 7.0, 6.0, 5.0, 4.0}, {4.0, 1.0, 3.5, 3.0, 2.5, 6.5, 5.0, 6.0, 3.5},
                                                {5.0, 2.0, 5.0, 3.0, 3.0, 6.0, 4.0, 7.0, 3.0}, {5.0, 2.0, 5.0, 3.0, 3.0, 6.0, 4.0, 7.0, 3.0}};

    this->create_mesh("mesh/cube_2x1.msh");
    this->read_input(eq_data_input);

    for (unsigned int i=0; i<n_times; i++) {  // time loop
        eq_data_->reallocate_cache();
        SingleValRef<double> ref_scalar(expected_scalar[i]);
        SingleValRef<arma::vec3> ref_vector(expected_vector[i]);
        SingleValRef<arma::mat33> ref_tensor(expected_tensor[i]);
        EXPECT_TRUE( eval_bulk_field(eq_data_->scalar_field, ref_scalar) );
        EXPECT_TRUE( eval_bulk_field(eq_data_->vector_field, ref_vector) );
        EXPECT_TRUE( eval_bulk_field(eq_data_->tensor_field, ref_tensor) );
        eq_data_->tg_.next_time();
    }
}
