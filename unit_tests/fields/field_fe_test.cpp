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
#include "fields/eval_points.hh"
#include "fields/eval_subset.hh"
#include "fields/field_value_cache.hh"
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




//class FieldFETest : public testing::Test {
//public:
//    typedef FieldFE<3, FieldValue<3>::Scalar > ScalarField;
//    typedef FieldFE<3, FieldValue<3>::VectorFixed > VecField;
//
//    virtual void SetUp() {
//    	this->mesh = nullptr;
//        // setup FilePath directories
//        FilePath::set_io_dirs(".",UNIT_TESTS_SRC_DIR,"",".");
//
//        Profiler::instance();
//        PetscInitialize(0,PETSC_NULL,PETSC_NULL,PETSC_NULL);
//    }
//
//    virtual void TearDown() {
//    	dh.reset();
//    	if (mesh != nullptr) delete mesh;
//    }
//
//    void create_mesh(std::string mesh_file_str) {
//        mesh = mesh_full_constructor("{ mesh_file=\"" + mesh_file_str + "\", optimize_mesh=false }");
//    }
//
//    void create_dof_handler(double val1, double val2, double val3) {
//        dh = std::make_shared<DOFHandlerMultiDim>(*mesh);
//        v.resize(3);
//        v.set(0, val1);
//        v.set(1, val2);
//        v.set(2, val3);
//        dof_values[0] = val1;
//        dof_values[1] = val2;
//        dof_values[2] = val3;
//    }
//
//    const FieldAlgoBaseInitData& init_data(std::string field_name) {
//    	static const FieldAlgoBaseInitData init_data(field_name, 0, UnitSI::dimensionless());
//    	return init_data;
//    }
//
//    static Input::Type::Record &get_input_type() {
//        return Input::Type::Record("Test","")
//            .declare_key("scalar", ScalarField::get_input_type(), Input::Type::Default::obligatory(),"" )
//            .declare_key("native_data", ScalarField::get_input_type(), Input::Type::Default::obligatory(),"" )
//            .close();
//    }
//
//    Mesh *mesh;
//    std::shared_ptr<DOFHandlerMultiDim> dh;
//    double dof_values[3];
//    VectorMPI v;
//
//};
//
//
//TEST_F(FieldFETest, scalar) {  // moved to FieldEvalFETest, set_fe_data_scalar
//    create_mesh("fields/one_element_2d.msh");
//    create_dof_handler(1, 2, 3);
//
//	MixedPtr<FE_P_disc> fe(1);
//    std::shared_ptr<DiscreteSpace> ds = std::make_shared<EqualOrderDiscreteSpace>(mesh, fe);
//    ScalarField field;
//
//    dh->distribute_dofs(ds);
//    field.set_fe_data(dh, v);
//    field.set_time(0.0);
//
//    Armor::array pts(3, 1);
//    pts.reinit(3);
//    pts.append(Armor::vec<3>({ 1, 1, 5 }));
//    pts.append(Armor::vec<3>({ 4, 0, 5 }));
//    pts.append(Armor::vec<3>({ 2, 3, 5 }));
//    vector<double> values(3);
//
//    // test values at vertices of the triangle
//    field.value_list( pts, mesh->element_accessor(0), values );
//    EXPECT_DOUBLE_EQ( dof_values[0], values[0] );
//    EXPECT_DOUBLE_EQ( dof_values[1], values[1] );
//    EXPECT_DOUBLE_EQ( dof_values[2], values[2] );
//
//    // test value at barycenter
//    EXPECT_DOUBLE_EQ( (dof_values[0]+dof_values[1]+dof_values[2])/3, field.value({ 7./3, 4./3, 5 }, mesh->element_accessor(0)) );
//}
//
//
//TEST_F(FieldFETest, vector) {  // moved to FieldEvalFETest, set_fe_data_vector
//    create_mesh("fields/one_element_2d.msh");
//    create_dof_handler(0, 0, 1);
//
//	MixedPtr<FE_RT0> fe;
//    std::shared_ptr<DiscreteSpace> ds = std::make_shared<EqualOrderDiscreteSpace>(mesh, fe);
//    VecField field;
//
//    dh->distribute_dofs(ds);
//    field.set_fe_data(dh, v);
//    field.set_time(0.0);
//
//    // The Raviart-Thomas function given by the following dofs
//    // is 3/7*(x-7/3, y-4/3, 0).
//
//    arma::vec3 result = { 2./7, 1./14, 0 };
//
//    EXPECT_NEAR( 0, arma::norm(result - field.value({ 3, 1.5, 5 }, mesh->element_accessor(0)), 2), 1e-15 );
//}
//
//
//string input = R"INPUT(
//{
//   scalar={
//       TYPE="FieldFE",
//       mesh_data_file="fields/simplest_cube_data.msh",
//       field_name="scalar"
//   }
//   native_data={
//       TYPE="FieldFE",
//       mesh_data_file="output/test_output_vtk_ascii_ref.vtu",
//       field_name="flow_data"
//   }
//}
//)INPUT";
//
//
//
//TEST_F(FieldFETest, scalar_from_input) { // equivalent with TEST(FieldFENewTest, scalar) - moved to FieldEvalFETest, input_msh
//    create_mesh("fields/simplest_cube_data.msh");
//
//    Input::ReaderToStorage reader( input, FieldFETest::get_input_type(), Input::FileFormat::format_JSON );
//    Input::Record rec=reader.get_root_interface<Input::Record>();
//
//    ScalarField field;
//    field.init_from_input(rec.val<Input::Record>("scalar"), init_data("scalar"));
//    field.set_mesh(mesh,false);
//    field.set_time(0.0);
//
//    Space<3>::Point point;
//    for(unsigned int i=0; i < mesh->n_elements(); i++) {
//        EXPECT_DOUBLE_EQ( (i+1)*0.1 , field.value(point, mesh->element_accessor(i)) );
//    }
//}
//
//
//TEST_F(FieldFETest, native_data) { // moved to FieldEvalFETest, native_data
//    create_mesh("fields/simplest_cube_3d.msh");
//
//    Input::ReaderToStorage reader( input, FieldFETest::get_input_type(), Input::FileFormat::format_JSON );
//    Input::Record rec=reader.get_root_interface<Input::Record>();
//
//    ScalarField field;
//    field.init_from_input(rec.val<Input::Record>("native_data"), init_data("native_data"));
//    field.set_mesh(mesh,false);
//    field.set_time(0.0);
//
//    Space<3>::Point point;
//    for(unsigned int i=0; i < mesh->n_elements(); i++) {
//        EXPECT_DOUBLE_EQ( i*0.2 , field.value(point, mesh->element_accessor(i)) );
//    }
//}


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
          default_value: 0.0
        vector_field: !FieldFE
          mesh_data_file: fields/simplest_cube_data.msh
          field_name: vector_fixed
          default_value: 0.0
        tensor_field: !FieldFE
          mesh_data_file: fields/simplest_cube_data.msh
          field_name: tensor_fixed
          default_value: 0.0
        enum_field: !FieldFE
          mesh_data_file: fields/simplest_cube_data.msh
          field_name: enum
          default_value: 0
        bc_scalar_field: !FieldFE
          mesh_data_file: fields/simplest_cube_data.msh
          field_name: scalar
          default_value: 0.0
        bc_vector_field: !FieldFE
          mesh_data_file: fields/simplest_cube_data.msh
          field_name: vector_fixed
          default_value: 0.0
        bc_tensor_field: !FieldFE
          mesh_data_file: fields/simplest_cube_data.msh
          field_name: tensor_fixed
          default_value: 0.0
        bc_enum_field: !FieldFE
          mesh_data_file: fields/simplest_cube_data.msh
          field_name: enum
          default_value: 0
        scalar_ref: !FieldFormula
          value: x+2*y+t
        vector_ref: !FieldFormula
          value: [x+2*y, y+2*z+0.5*t, z+2*x+t]
        tensor_ref: !FieldFormula
          value: [2*x+y, 2*y+z+0.5*t, 2*z+x+t]
        bc_scalar_ref: !FieldFormula
          value: x+y+z+t
        bc_vector_ref: !FieldFormula
          value: [3*x, 3*y+0.5*t, 3*z+t]
        bc_tensor_ref: !FieldFormula
          value: [x+y+z, 2*(x+y+z)+0.5*t, 3*(x+y+z)+t]
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
        EXPECT_TRUE( eval_boundary_field(eq_data_->bc_scalar_field, eq_data_->bc_scalar_ref, 3, 0) );
        EXPECT_TRUE( eval_boundary_field(eq_data_->bc_vector_field, eq_data_->bc_vector_ref, 3, 0) );
        EXPECT_TRUE( eval_boundary_field(eq_data_->bc_tensor_field, eq_data_->bc_tensor_ref, 3, 0) );
        EXPECT_TRUE( eval_boundary_field(eq_data_->bc_enum_field, j+1, 3, 0) );
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
          value: x+2*y
        vector_ref: !FieldFormula
          value: [x+2*y, y+2*z, z+2*x]
        tensor_ref: !FieldFormula
          value: [2*x+y, 2*y+z, 2*z+x]
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
          default_value: 0.0
          read_time_shift: 1.0
        scalar_ref: !FieldFormula
          value: x+2*y+(t+1)
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
      - region: ALL
        time: 0.0
        scalar_field: !FieldFE
          mesh_data_file: fields/simplest_cube_data.msh
          field_name: scalar
          unit: "const; const=0.1*m"
          default_value: 0.0
        bc_scalar_field: !FieldFE
          mesh_data_file: fields/simplest_cube_data.msh
          field_name: scalar
          unit: "const; const=0.1*m"
          default_value: 0.0
        scalar_ref: !FieldFormula
          value: 0.1*(x+2*y+t)
        bc_scalar_ref: !FieldFormula
          value: 0.1*(x+y+z+t)
    )YAML";

    this->create_mesh("mesh/simplest_cube.msh");
    this->read_input(eq_data_input);

    for (unsigned int j=0; j<2; j++) {  // time loop
        eq_data_->reallocate_cache();

        // BULK field
    	FieldRef<ScalarField> ref_scalar(eq_data_->scalar_ref);
        EXPECT_TRUE( eval_bulk_field(eq_data_->scalar_field, ref_scalar) );

        // BOUNDARY field
        EXPECT_TRUE( eval_boundary_field(eq_data_->bc_scalar_field, eq_data_->bc_scalar_ref, 3, 0) );
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
          default_value: 0.0
          interpolation: identic_mesh
        vector_field: !FieldFE
          mesh_data_file: fields/identic_mesh_data.msh
          field_name: vector_fixed
          default_value: 0.0
          interpolation: identic_mesh
        bc_scalar_field: !FieldFE
          mesh_data_file: fields/identic_mesh_data.msh
          field_name: scalar
          default_value: 0.0
          interpolation: identic_mesh
        bc_vector_field: !FieldFE
          mesh_data_file: fields/identic_mesh_data.msh
          field_name: vector_fixed
          default_value: 0.0
          interpolation: identic_mesh
        scalar_ref: !FieldFormula
          value: x+2*y+t
        vector_ref: !FieldFormula
          value: [x+2*y, y+2*z+0.5*t, z+2*x+t]
        bc_scalar_ref: !FieldFormula
          value: x+y+z+t
        bc_vector_ref: !FieldFormula
          value: [3*x, 3*y+0.5*t, 3*z+t]
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
        EXPECT_TRUE( eval_boundary_field(eq_data_->bc_scalar_field, eq_data_->bc_scalar_ref, 3, 0) );
        EXPECT_TRUE( eval_boundary_field(eq_data_->bc_vector_field, eq_data_->bc_vector_ref, 3, 0) );
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
        bc_vector_field: !FieldFE
          mesh_data_file: fields/interpolation_rectangle.msh
          field_name: vector_fixed
          default_value: 0.0
          #interpolation: P0_intersection
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
