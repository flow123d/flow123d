/*
 * field_python_test.cpp
 *
 *  Created on: Aug 30, 2012
 *      Author: jb
 */


#define TEST_USE_MPI
#define FEAL_OVERRIDE_ASSERTS
#include <flow_gtest_mpi.hh>
#include <mesh_constructor.hh>
#include "arma_expect.hh"

#include "field_eval_base.hh"              // for FieldEvalBaseTest
#include "fields/eval_points.hh"
#include "fields/eval_subset.hh"
#include "fields/field_value_cache.hh"
#include "fields/field_values.hh"
#include "fields/field_set.hh"
#include "tools/unit_si.hh"
#include "fields/bc_field.hh"
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
#include "system/python_loader.hh"


class FieldEvalPythonTest : public FieldEvalBaseTest {
public:

    class EqData : public FieldEvalBaseTest::EqData {
        EqData() : FieldEvalBaseTest::EqData() {}
    };

    FieldEvalPythonTest() : FieldEvalBaseTest() {}

    ~FieldEvalPythonTest() {}
};


TEST_F(FieldEvalPythonTest, evaluate) {
    string eq_data_input = R"YAML(
    data:
      - region: 3D left
        time: 0.0
        scalar_field: !FieldPython
          source_file: ./fields/field_python_test.py
          class: FieldPythonTest1
          used_fields: ['X']
        vector_field: !FieldPython
          source_file: ./fields/field_python_test.py
          class: FieldPythonTest1
          used_fields: ['X']
        tensor_field: !FieldPython
          source_file: ./fields/field_python_test.py
          class: FieldPythonTest1
          used_fields: ['X']
        scalar_ref: !FieldFormula
          value: X[0]
        vector_ref: !FieldFormula
          value: "[X[0], 2*X[0], 0.5]"
        tensor_ref: !FieldFormula
          value: "[ [X[0], 0.2, 0.3], [0.2, 0.4, 0.5], [0.3, 0.5, 0.6] ]"
      - region: 3D right
        time: 0.0
        scalar_field: !FieldPython
          source_file: ./fields/field_python_test.py
          class: FieldPythonTest2
          used_fields: ['X']
        vector_field:  !FieldPython
          source_file: ./fields/field_python_test.py
          class: FieldPythonTest2
          used_fields: ['X']
        tensor_field: !FieldPython
          source_file: ./fields/field_python_test.py
          class: FieldPythonTest2
          used_fields: ['X']
        scalar_ref: !FieldFormula
          value: X[1]
        vector_ref: !FieldFormula
          value: "[X[1], 2*X[1], 0.5]"
        tensor_ref: !FieldFormula
          value: "[ [X[1], 2.2, 2.3], [2.2, 2.4, 2.5], [2.3, 2.5, 2.6] ]"
    )YAML";

    this->create_mesh("mesh/cube_2x1.msh");
    this->read_input(eq_data_input);

    eq_data_->reallocate_cache();

    FieldRef<ScalarField> ref_scalar(eq_data_->scalar_ref);
    FieldRef<VectorField> ref_vector(eq_data_->vector_ref);
    FieldRef<TensorField> ref_tensor(eq_data_->tensor_ref);
    EXPECT_TRUE( eval_bulk_field(eq_data_->scalar_field, ref_scalar) );
    EXPECT_TRUE( eval_bulk_field(eq_data_->vector_field, ref_vector) );
    EXPECT_TRUE( eval_bulk_field(eq_data_->tensor_field, ref_tensor) );

}


TEST_F(FieldEvalPythonTest, exc_nonexist_file) {
    string eq_data_input = R"YAML(
    data:
      - region: BULK
        time: 0.0
        scalar_field: !FieldPython
          source_file: fields/field_python_asm.py
          class: SomeClass
          used_fields: ['X']
    )YAML";

    this->create_mesh("mesh/cube_2x1.msh");
    this->read_input(eq_data_input);
    EXPECT_THROW_WHAT( { eq_data_->reallocate_cache(); }, FilePath::ExcFileOpen,
        "Can not open file:");
}


TEST_F(FieldEvalPythonTest, exc_nonexist_class) {
    string eq_data_input = R"YAML(
    data:
      - region: BULK
        time: 0.0
        scalar_field: !FieldPython
          source_file: fields/field_python_test.py
          class: NonExistClass
          used_fields: ['X']
    )YAML";

    this->create_mesh("mesh/cube_2x1.msh");
    this->read_input(eq_data_input);
    EXPECT_THROW_WHAT( { eq_data_->reallocate_cache(); }, PythonLoader::ExcPythonError,
        "has no attribute 'NonExistClass'");
}


TEST_F(FieldEvalPythonTest, exc_nonexist_function) {
    string eq_data_input = R"YAML(
    data:
      - region: BULK
        time: 0.0
        scalar_field: !FieldPython
          source_file: fields/field_python_test.py
          class: FieldPythonTest3
          used_fields: ['X']
    )YAML";

    this->create_mesh("mesh/cube_2x1.msh");
    this->read_input(eq_data_input);
    eq_data_->reallocate_cache();

    SingleValRef<double> ref_scalar(0.0);
    EXPECT_THROW_WHAT( { eval_bulk_field(eq_data_->scalar_field, ref_scalar); }, PythonLoader::ExcPythonError,
        "'NoneType' object is not subscriptable");
}

