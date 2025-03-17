/*
 * field_fe_test.cpp
 *
 *  Created on: Jan 25, 2013
 *      Author: jb
 */



#define TEST_USE_MPI
#define FEAL_OVERRIDE_ASSERTS
#include <flow_gtest_mpi.hh>
#include <mesh_constructor.hh>
#include <limits>
#include "arma_expect.hh"


#include "field_eval_base.hh"
#include "fem/eval_points.hh"
#include "fem/eval_subset.hh"
#include "fem/element_cache_map.hh"
#include "fields/field_values.hh"
#include "fields/field_set.hh"
#include "fields/field_fe.hh"
#include "tools/unit_si.hh"
#include "fields/bc_field.hh"
#include "quadrature/quadrature.hh"
#include "quadrature/quadrature_lib.hh"
#include "fem/dofhandler.hh"
#include "fem/dh_cell_accessor.hh"
#include "fem/fe_p.hh"
#include "fem/fe_rt.hh"
#include "mesh/mesh.h"
#include "mesh/accessors.hh"
#include "input/input_type.hh"
#include "input/accessors.hh"
#include "input/reader_to_storage.hh"
#include "system/sys_profiler.hh"
#include "la/vector_mpi.hh"
#include "tools/mixed.hh"




class FieldEvalFETest : public FieldEvalBaseTest {
public:

    class EqData : public FieldEvalBaseTest::EqData {
        EqData() : FieldEvalBaseTest::EqData() {}
    };

    FieldEvalFETest() : FieldEvalBaseTest() {}

    ~FieldEvalFETest() {}

    /// Specialization of set_fe_data unit tests
    void set_dof_values(std::vector<double> vals) {
        v.resize(vals.size());
        dof_values.resize(vals.size());
        for (unsigned int i=0; i<vals.size(); ++i) {
            v.set(i, vals[i]);
            dof_values[i] = vals[i];
        }

        eq_data_->set_mesh(*mesh_);
    }


    std::vector<double> dof_values;           ///< used in test set_fe_data
    VectorMPI v;                              ///< used in test set_fe_data
};





TEST_F(FieldEvalFETest, input_msh) {
    string eq_data_input = R"YAML(
    data:
      - region: ALL
        time: 0.0
        scalar_field: !FieldFE
          mesh_data_file: fields/simplest_cube_data.msh
          field_name: scalar
        vector_field: !FieldFE
          mesh_data_file: fields/simplest_cube_data.msh
          field_name: vector_fixed
        tensor_field: !FieldFE
          mesh_data_file: fields/simplest_cube_data.msh
          field_name: tensor_fixed
        enum_field: !FieldFE
          mesh_data_file: fields/simplest_cube_data.msh
          field_name: enum
        scalar_ref: !FieldFormula
          value: X[0]+2*X[1]+t
        vector_ref: !FieldFormula
          value: "[X[0]+2*X[1], X[1]+2*X[2]+0.5*t, X[2]+2*X[0]+t]"
        tensor_ref: !FieldFormula
          value: "[ [2*X[0]+X[1], 0, 0], [0, 2*X[1]+X[2]+0.5*t, 0], [0, 0, 2*X[2]+X[0]+t] ]"
      - region: [".top side", ".bottom side"]
        time: 0.0
        bc_scalar_field: !FieldFE
          mesh_data_file: fields/simplest_cube_data.msh
          field_name: scalar
          is_boundary: true
        bc_vector_field: !FieldFE
          mesh_data_file: fields/simplest_cube_data.msh
          field_name: vector_fixed
          is_boundary: true
        bc_tensor_field: !FieldFE
          mesh_data_file: fields/simplest_cube_data.msh
          field_name: tensor_fixed
          is_boundary: true
        bc_enum_field: !FieldFE
          mesh_data_file: fields/simplest_cube_data.msh
          field_name: enum
          is_boundary: true
        bc_scalar_ref: !FieldFormula
          value: X[0]+X[1]+X[2]+t
        bc_vector_ref: !FieldFormula
          value: "[3*X[0], 3*X[1]+0.5*t, 3*X[2]+t]"
        bc_tensor_ref: !FieldFormula
          value: "[ [X[0]+X[1]+X[2], 0, 0], [0, 2*(X[0]+X[1]+X[2])+0.5*t, 0], [0, 0, 3*(X[0]+X[1]+X[2])+t] ]"
    )YAML";

    this->create_mesh("mesh/simplest_cube.msh");
    this->read_input(eq_data_input);

    for (unsigned int j=0; j<2; j++) {  // time loop
    	eq_data_->reallocate_cache();

        // BULK fields
    	FieldRef<ScalarField> ref_scalar(eq_data_->scalar_ref);
    	FieldRef<VectorField> ref_vector(eq_data_->vector_ref);
    	FieldRef<TensorField> ref_tensor(eq_data_->tensor_ref);
    	SingleValRef<unsigned int> ref_enum(j);
        EXPECT_TRUE( eval_bulk_field(eq_data_->scalar_field, ref_scalar) );
        EXPECT_TRUE( eval_bulk_field(eq_data_->vector_field, ref_vector) );
        EXPECT_TRUE( eval_bulk_field(eq_data_->tensor_field, ref_tensor) );
        EXPECT_TRUE( eval_bulk_field(eq_data_->enum_field, ref_enum) );

        // BOUNDARY fields
        FieldRef<ScalarField> ref_bc_scalar(eq_data_->bc_scalar_ref);
        FieldRef<VectorField> ref_bc_vector(eq_data_->bc_vector_ref);
        FieldRef<TensorField> ref_bc_tensor(eq_data_->bc_tensor_ref);
        SingleValRef<unsigned int> ref_bc_enum(j+1);
        EXPECT_TRUE( eval_boundary_field(eq_data_->bc_scalar_field, ref_bc_scalar, 3, 0) );
        EXPECT_TRUE( eval_boundary_field(eq_data_->bc_vector_field, ref_bc_vector, 3, 0) );
        EXPECT_TRUE( eval_boundary_field(eq_data_->bc_tensor_field, ref_bc_tensor, 3, 0) );
        EXPECT_TRUE( eval_boundary_field(eq_data_->bc_enum_field, ref_bc_enum, 3, 0) );
        eq_data_->tg_.next_time();
    }
}


TEST_F(FieldEvalFETest, input_vtk) {
    string eq_data_input = R"YAML(
    data:
      - region: ALL
        time: 0.0
        scalar_field: !FieldFE
          mesh_data_file: fields/vtk_ascii_data.vtu
          field_name: scalar_field
        vector_field: !FieldFE
          mesh_data_file: fields/vtk_ascii_data.vtu
          field_name: vector_field
        tensor_field: !FieldFE
          mesh_data_file: fields/vtk_ascii_data.vtu
          field_name: tensor_field
        scalar_ref: !FieldFormula
          value: X[0]+2*X[1]
        vector_ref: !FieldFormula
          value: "[X[0]+2*X[1], X[1]+2*X[2], X[2]+2*X[0]]"
        tensor_ref: !FieldFormula
          value: "[ [2*X[0]+X[1], 0, 0], [0, 2*X[1]+X[2], 0], [0, 0, 2*X[2]+X[0]] ]"
    )YAML";

    this->create_mesh("mesh/simplest_cube.msh");
    this->read_input(eq_data_input);

    eq_data_->reallocate_cache();
	FieldRef<ScalarField> ref_scalar(eq_data_->scalar_ref);
	FieldRef<VectorField> ref_vector(eq_data_->vector_ref);
	FieldRef<TensorField> ref_tensor(eq_data_->tensor_ref);
    EXPECT_TRUE( eval_bulk_field(eq_data_->scalar_field, ref_scalar) );
    EXPECT_TRUE( eval_bulk_field(eq_data_->vector_field, ref_vector) );
    EXPECT_TRUE( eval_bulk_field(eq_data_->tensor_field, ref_tensor) );
}


TEST_F(FieldEvalFETest, time_shift) {
    string eq_data_input = R"YAML(
    data:
      - region: ALL
        time: 0.0
        scalar_field: !FieldFE
          mesh_data_file: fields/simplest_cube_data.msh
          field_name: scalar
          read_time_shift: 1.0
        scalar_ref: !FieldFormula
          value: X[0]+2*X[1]+(t+1)
    )YAML";

    this->create_mesh("mesh/simplest_cube.msh");
    this->read_input(eq_data_input);

    for (unsigned int j=0; j<2; j++) {  // time loop
        eq_data_->reallocate_cache();
    	FieldRef<ScalarField> ref_scalar(eq_data_->scalar_ref);
        EXPECT_TRUE( eval_bulk_field(eq_data_->scalar_field, ref_scalar) );
        eq_data_->tg_.next_time();
    }
}


TEST_F(FieldEvalFETest, default_values) {
    string eq_data_input = R"YAML(
    data:
      - region: ALL
        time: 0.0
        vector_field: !FieldFE
          mesh_data_file: fields/simplest_cube_data.msh
          field_name: unset_bulk
          default_value: 0.1
    )YAML";

    arma::vec3 expected_vector = arma::vec3("0.1 0.1 0.1");

    this->create_mesh("mesh/simplest_cube.msh");
    this->read_input(eq_data_input);

    for (unsigned int j=0; j<2; j++) {  // time loop
        eq_data_->reallocate_cache();
    	SingleValRef<arma::vec3> ref_vector(expected_vector);
        EXPECT_TRUE( eval_bulk_field(eq_data_->vector_field, ref_vector) );
        eq_data_->tg_.next_time();
    }
}


TEST_F(FieldEvalFETest, unit_conversion) {
    string eq_data_input = R"YAML(
    data:
      - region: BULK
        time: 0.0
        scalar_field: !FieldFE
          mesh_data_file: fields/simplest_cube_data.msh
          field_name: scalar
          unit: "const; const=0.1*m"
        scalar_ref: !FieldFormula
          value: 0.1*(X[0]+2*X[1]+t)
      - region: [".top side", ".bottom side"]
        time: 0.0
        bc_scalar_field: !FieldFE
          mesh_data_file: fields/simplest_cube_data.msh
          field_name: scalar
          unit: "const; const=0.1*m"
          is_boundary: true
        bc_scalar_ref: !FieldFormula
          value: 0.1*(X[0]+X[1]+X[2]+t)
    )YAML";

    this->create_mesh("mesh/simplest_cube.msh");
    this->read_input(eq_data_input);

    for (unsigned int j=0; j<2; j++) {  // time loop
        eq_data_->reallocate_cache();

        // BULK field
    	FieldRef<ScalarField> ref_scalar(eq_data_->scalar_ref);
        EXPECT_TRUE( eval_bulk_field(eq_data_->scalar_field, ref_scalar) );

        // BOUNDARY field
    	FieldRef<ScalarField> ref_bc_scalar(eq_data_->bc_scalar_ref);
        EXPECT_TRUE( eval_boundary_field(eq_data_->bc_scalar_field, ref_bc_scalar, 3, 0) );
        eq_data_->tg_.next_time();
    }
}


TEST_F(FieldEvalFETest, identic_mesh) {
    string eq_data_input = R"YAML(
    data:
      - region: ALL
        time: 0.0
        scalar_field: !FieldFE
          mesh_data_file: fields/identic_mesh_data.msh
          field_name: scalar
          interpolation: identic_mesh
        vector_field: !FieldFE
          mesh_data_file: fields/identic_mesh_data.msh
          field_name: vector_fixed
          interpolation: identic_mesh
        scalar_ref: !FieldFormula
          value: X[0]+2*X[1]+t
        vector_ref: !FieldFormula
          value: '[X[0]+2*X[1], X[1]+2*X[2]+0.5*t, X[2]+2*X[0]+t]'
      - region: [".top side", ".bottom side"]
        time: 0.0
        bc_scalar_field: !FieldFE
          mesh_data_file: fields/identic_mesh_data.msh
          field_name: scalar
          interpolation: identic_mesh
          is_boundary: true
        bc_vector_field: !FieldFE
          mesh_data_file: fields/identic_mesh_data.msh
          field_name: vector_fixed
          interpolation: identic_mesh
          is_boundary: true
        bc_scalar_ref: !FieldFormula
          value: X[0]+X[1]+X[2]+t
        bc_vector_ref: !FieldFormula
          value: '[3*X[0], 3*X[1]+0.5*t, 3*X[2]+t]'
    )YAML";

    this->create_mesh("mesh/simplest_cube.msh");
    this->read_input(eq_data_input);

    for (unsigned int j=0; j<2; j++) {  // time loop
        eq_data_->reallocate_cache();

        // BULK field
    	FieldRef<ScalarField> ref_scalar(eq_data_->scalar_ref);
    	FieldRef<VectorField> ref_vector(eq_data_->vector_ref);
        EXPECT_TRUE( eval_bulk_field(eq_data_->scalar_field, ref_scalar) );
        EXPECT_TRUE( eval_bulk_field(eq_data_->vector_field, ref_vector) );

        // BOUNDARY field
    	FieldRef<ScalarField> ref_bc_scalar(eq_data_->bc_scalar_ref);
    	FieldRef<VectorField> ref_bc_vector(eq_data_->bc_vector_ref);
        EXPECT_TRUE( eval_boundary_field(eq_data_->bc_scalar_field, ref_bc_scalar, 3, 0) );
        EXPECT_TRUE( eval_boundary_field(eq_data_->bc_vector_field, ref_bc_vector, 3, 0) );
        eq_data_->tg_.next_time();
    }
}


TEST_F(FieldEvalFETest, interpolation_gauss) {
    string eq_data_input = R"YAML(
    data:
      - region: ALL
        time: 0.0
        scalar_field: !FieldFE
          mesh_data_file: fields/interpolation_rectangle.msh
          field_name: scalar
          default_value: 0.0
          #interpolation: P0_gauss
        vector_field: !FieldFE
          mesh_data_file: fields/interpolation_rectangle.msh
          field_name: vector_fixed
          default_value: 0.0
          #interpolation: P0_gauss
    )YAML";

    std::vector< std::vector<double> >     expected_scalars = { {0.25, 0.15, 0.25, 0.35}, {0.75, 0.65, 0.75, 0.85} };
    std::vector< std::vector<arma::vec3> > expected_vectors = { {"2.5 3.5 4.5", "1.5 2.5 3.5", "2.5 3.5 4.5", "3.5 4.5 5.5"},
                                                                {"5.5 6.5 7.5", "4.5 5.5 6.5", "5.5 6.5 7.5", "6.5 7.5 8.5"} };

    this->create_mesh("fields/interpolation_rect_small.msh");
    this->read_input(eq_data_input);

    for (unsigned int j=0; j<2; j++) {  // time loop
        eq_data_->reallocate_cache();
        VecRef<double> ref_scalar(expected_scalars[j]);
        VecRef<arma::vec3> ref_vector(expected_vectors[j]);
        EXPECT_TRUE( eval_bulk_field(eq_data_->scalar_field, ref_scalar) );
        EXPECT_TRUE( eval_bulk_field(eq_data_->vector_field, ref_vector) );
        eq_data_->tg_.next_time();
    }
}


TEST_F(FieldEvalFETest, interpolation_gauss_unit_conversion) {
    string eq_data_input = R"YAML(
    data:
      - region: ALL
        time: 0.0
        scalar_field: !FieldFE
          mesh_data_file: fields/interpolation_rectangle.msh
          field_name: scalar
          default_value: 0.0
          unit: "conv; conv=10*m"
          #interpolation: P0_gauss
    )YAML";

    std::vector< std::vector<double> > expected_scalars = { {2.5, 1.5, 2.5, 3.5}, {7.5, 6.5, 7.5, 8.5} };

    this->create_mesh("fields/interpolation_rect_small.msh");
    this->read_input(eq_data_input);

    for (unsigned int j=0; j<2; j++) {  // time loop
        eq_data_->reallocate_cache();
        VecRef<double> ref_scalar(expected_scalars[j]);
        EXPECT_TRUE( eval_bulk_field(eq_data_->scalar_field, ref_scalar) );
        eq_data_->tg_.next_time();
    }
}


TEST_F(FieldEvalFETest, interpolation_1d_2d) { // TODO fix bdr
    string eq_data_input = R"YAML(
    data:
      - region: ALL
        time: 0.0
        bc_scalar_field: !FieldFE
          mesh_data_file: fields/interpolation_rectangle.msh
          field_name: scalar
          default_value: 0.0
          #interpolation: P0_intersection
          is_boundary: true
        bc_vector_field: !FieldFE
          mesh_data_file: fields/interpolation_rectangle.msh
          field_name: vector_fixed
          default_value: 0.0
          #interpolation: P0_intersection
          is_boundary: true
    )YAML";

    std::vector< std::vector<double> >     expected_scalars = { {0.25, 0.15, 0.25, 0.35}, {0.75, 0.65, 0.75, 0.85} };
    std::vector< std::vector<arma::vec3> > expected_vectors = {{"2.5 3.5 4.5", "1.5 2.5 3.5", "2.5 3.5 4.5", "3.5 4.5 5.5"},
                                                               {"5.5 6.5 7.5", "4.5 5.5 6.5", "5.5 6.5 7.5", "6.5 7.5 8.5"} };

    this->create_mesh("fields/interpolation_rect_small.msh");
    this->read_input(eq_data_input);

    for (unsigned int j=0; j<2; j++) {
        eq_data_->reallocate_cache();

        for(unsigned int i=0; i < mesh_->n_elements(); i++) {  // time loop
            SingleValRef<double> ref_scalar(expected_scalars[j][i]);
            SingleValRef<arma::vec3> ref_vector(expected_vectors[j][i]);
            EXPECT_TRUE( eval_boundary_field(eq_data_->bc_scalar_field, ref_scalar, i, i) );
            EXPECT_TRUE( eval_boundary_field(eq_data_->bc_vector_field, ref_vector, i, i) );
        }
        eq_data_->tg_.next_time();
    }
}


TEST_F(FieldEvalFETest, native) {
    string eq_data_input = R"YAML(
    data:
      - region: ALL
        time: 0.0
        scalar_field: !FieldFE
          mesh_data_file: output/test_output_vtk_ascii_ref.vtu
          field_name: flow_data
    )YAML";

    std::vector<double> expected_scalars = { 0.0, 0.2, 0.4, 0.6, 0.8, 1.0 };

    this->create_mesh("fields/simplest_cube_3d.msh");
    this->read_input(eq_data_input);

    eq_data_->reallocate_cache();
    VecRef<double> ref_scalar(expected_scalars);
    EXPECT_TRUE( eval_bulk_field(eq_data_->scalar_field, ref_scalar) );
}


TEST_F(FieldEvalFETest, set_fe_data_scalar) {
    typedef FieldFE<3, FieldValue<3>::Scalar > ScalarFieldFE;
    this->create_mesh("fields/one_element_2d.msh");
    this->set_dof_values( {0.5} );

    MixedPtr<FE_P_disc> fe(0);
    std::shared_ptr<DiscreteSpace> ds = std::make_shared<EqualOrderDiscreteSpace>(mesh_, fe);
    dh_->distribute_dofs(ds);

    std::shared_ptr<ScalarFieldFE> fe_field = std::make_shared<ScalarFieldFE>();
    fe_field->set_fe_data(dh_, v);
    eq_data_->scalar_field.set(fe_field, 0.0);

    eq_data_->reallocate_cache();
	SingleValRef<double> ref_scalar(0.5);
    EXPECT_TRUE( eval_bulk_field(eq_data_->scalar_field, ref_scalar) );
}


TEST_F(FieldEvalFETest, set_fe_data_vector) {
    typedef FieldFE<3, FieldValue<3>::VectorFixed > VectorFieldFE;
    this->create_mesh("fields/one_element_2d.msh");
    this->set_dof_values( {0.5, 1.5, 2.5} );

    MixedPtr<FE_RT0> fe;
    std::shared_ptr<DiscreteSpace> ds = std::make_shared<EqualOrderDiscreteSpace>(mesh_, fe);
    dh_->distribute_dofs(ds);

    std::shared_ptr<VectorFieldFE> fe_field = std::make_shared<VectorFieldFE>();
    fe_field->set_fe_data(dh_, v);
    eq_data_->vector_field.set(fe_field, 0.0);

    eq_data_->reallocate_cache();
    arma::vec3 expected = { 1./7, 2./7, 0.0 };
    SingleValRef<arma::vec3> ref_vector(expected);
    EXPECT_TRUE( eval_bulk_field(eq_data_->vector_field, ref_vector) );
}
