/*
 * field_elementwise_test.cpp
 *
 *  Created on: Jan 25, 2013
 *      Author: jb
 */



#define TEST_USE_MPI
#define FEAL_OVERRIDE_ASSERTS
#include <flow_gtest_mpi.hh>
#include <mesh_constructor.hh>


#include "fields/field_fe.hh"
#include "la/vector_mpi.hh"
#include "fields/fe_value_handler.hh"
#include "tools/unit_si.hh"
#include "input/input_type.hh"
#include "input/accessors.hh"

#include "fem/dofhandler.hh"
#include "fem/fe_p.hh"
#include "quadrature/quadrature.hh"
#include "fem/mapping_p1.hh"
#include "fem/fe_values.hh"
#include "fem/fe_rt.hh"

#include "system/sys_profiler.hh"

#include "mesh/mesh.h"
#include "io/msh_gmshreader.h"
#include "io/reader_cache.hh"
#include "tools/mixed.hh"




class FieldFETest : public testing::Test {
public:
    typedef FieldFE<3, FieldValue<3>::Scalar > ScalarField;
    typedef FieldFE<3, FieldValue<3>::VectorFixed > VecField;

    virtual void SetUp() {
    	this->mesh = nullptr;
        // setup FilePath directories
        FilePath::set_io_dirs(".",UNIT_TESTS_SRC_DIR,"",".");

        Profiler::instance();
        PetscInitialize(0,PETSC_NULL,PETSC_NULL,PETSC_NULL);
    }

    virtual void TearDown() {
    	dh.reset();
    	if (mesh != nullptr) delete mesh;
    }

    void create_mesh(std::string mesh_file_str) {
        mesh = mesh_full_constructor("{mesh_file=\"" + mesh_file_str + "\"}");
    }

    void create_dof_handler(double val1, double val2, double val3) {
        dh = std::make_shared<DOFHandlerMultiDim>(*mesh);
        v.resize(3);
        v[0] = val1;
        v[1] = val2;
        v[2] = val3;
        dof_values[0] = val1;
        dof_values[1] = val2;
        dof_values[2] = val3;
    }

    const FieldAlgoBaseInitData& init_data(std::string field_name) {
    	static const FieldAlgoBaseInitData init_data(field_name, 0, UnitSI::dimensionless());
    	return init_data;
    }

    static Input::Type::Record &get_input_type() {
        return Input::Type::Record("Test","")
            .declare_key("scalar", ScalarField::get_input_type(), Input::Type::Default::obligatory(),"" )
            .declare_key("native_data", ScalarField::get_input_type(), Input::Type::Default::obligatory(),"" )
            .close();
    }

    Mesh *mesh;
    std::shared_ptr<DOFHandlerMultiDim> dh;
    double dof_values[3];
    VectorMPI v;

};


// TODO Fix these tests after improving DOF handler
TEST_F(FieldFETest, scalar) {
    create_mesh("fields/one_element_2d.msh");
    create_dof_handler(1, 2, 3);

	MixedPtr<FE_P_disc> fe(1);
    std::shared_ptr<DiscreteSpace> ds = std::make_shared<EqualOrderDiscreteSpace>(mesh, fe);
    ScalarField field;

    dh->distribute_dofs(ds);
    field.set_fe_data(dh, 0, v);
    field.set_time(0.0);

    Armor::array pts(3, 1);
    pts.reinit(3);
    pts.append(Armor::vec<3>({ 1, 1, 5 }));
    pts.append(Armor::vec<3>({ 4, 0, 5 }));
    pts.append(Armor::vec<3>({ 2, 3, 5 }));
    vector<double> values(3);

    // test values at vertices of the triangle
    field.value_list( pts, mesh->element_accessor(0), values );
    EXPECT_DOUBLE_EQ( dof_values[0], values[0] );
    EXPECT_DOUBLE_EQ( dof_values[1], values[1] );
    EXPECT_DOUBLE_EQ( dof_values[2], values[2] );

    // test value at barycenter
    EXPECT_DOUBLE_EQ( (dof_values[0]+dof_values[1]+dof_values[2])/3, field.value({ 7./3, 4./3, 5 }, mesh->element_accessor(0)) );
}


TEST_F(FieldFETest, vector) {
    create_mesh("fields/one_element_2d.msh");
    create_dof_handler(0, 0, 1);

	MixedPtr<FE_RT0> fe;
    std::shared_ptr<DiscreteSpace> ds = std::make_shared<EqualOrderDiscreteSpace>(mesh, fe);
    VecField field;

    dh->distribute_dofs(ds);
    field.set_fe_data(dh, 0, v);
    field.set_time(0.0);

    // The Raviart-Thomas function given by the following dofs
    // is 3/7*(x-7/3, y-4/3, 0).

    arma::vec3 result = { 2./7, 1./14, 0 };

    EXPECT_NEAR( 0, arma::norm(result - field.value({ 3, 1.5, 5 }, mesh->element_accessor(0)), 2), 1e-15 );
}


string input = R"INPUT(
{   
   scalar={
       TYPE="FieldFE",
       mesh_data_file="fields/simplest_cube_data.msh",
       field_name="scalar"
   }
   native_data={
       TYPE="FieldFE",
       mesh_data_file="output/test_output_vtk_ascii_ref.vtu",
       field_name="flow_data"
   }
}
)INPUT";



TEST_F(FieldFETest, scalar_from_input) {
    create_mesh("fields/simplest_cube_data.msh");

    Input::ReaderToStorage reader( input, FieldFETest::get_input_type(), Input::FileFormat::format_JSON );
    Input::Record rec=reader.get_root_interface<Input::Record>();

    ScalarField field;
    field.init_from_input(rec.val<Input::Record>("scalar"), init_data("scalar"));
    field.set_mesh(mesh,false);
    field.set_time(0.0);

    Space<3>::Point point;
    for(unsigned int i=0; i < mesh->n_elements(); i++) {
        EXPECT_DOUBLE_EQ( (i+1)*0.1 , field.value(point, mesh->element_accessor(i)) );
    }
}


TEST_F(FieldFETest, native_data) {
    create_mesh("fields/simplest_cube_3d.msh");

    Input::ReaderToStorage reader( input, FieldFETest::get_input_type(), Input::FileFormat::format_JSON );
    Input::Record rec=reader.get_root_interface<Input::Record>();

    ScalarField field;
    field.init_from_input(rec.val<Input::Record>("native_data"), init_data("native_data"));
    field.set_mesh(mesh,false);
    field.set_time(0.0);

    Space<3>::Point point;
    for(unsigned int i=0; i < mesh->n_elements(); i++) {
        EXPECT_DOUBLE_EQ( i*0.2 , field.value(point, mesh->element_accessor(i)) );
    }
}


/*******************************************************************************
 *                                                                             *
 *     New tests of Elementwise and Interpolation P0 replaced with FieldFE     *
 *                                                                             *
 *******************************************************************************/

string elem_input = R"YAML(
##### tests of elementwise
scalar: !FieldFE
  mesh_data_file: fields/simplest_cube_data.msh
  field_name: scalar
  default_value: 0.0
scalar_unit_conversion: !FieldFE
  mesh_data_file: fields/simplest_cube_data.msh
  field_name: scalar
  unit: "const; const=100*m^0"
  default_value: 0.0
scalar_time_shift: !FieldFE
  mesh_data_file: fields/simplest_cube_data.msh
  field_name: scalar
  default_value: 0.0
  read_time_shift: 1.0
enum: !FieldFE
  mesh_data_file: fields/simplest_cube_data.msh
  field_name: enum
  default_value: 0
vector_fixed: !FieldFE
  mesh_data_file: fields/simplest_cube_data.msh
  field_name: vector_fixed
  default_value: 0.0
tensor_fixed: !FieldFE
  mesh_data_file: fields/simplest_cube_data.msh
  field_name: tensor_fixed
  default_value: 0.0
vtk_scalar: !FieldFE
  mesh_data_file: fields/vtk_ascii_data.vtu
  field_name: scalar_field
vtk_vector: !FieldFE
  mesh_data_file: fields/vtk_binary_data.vtu
  field_name: vector_field
vtk_tensor: !FieldFE
  mesh_data_file: fields/vtk_compressed_data.vtu
  field_name: tensor_field
default_values: !FieldFE
  mesh_data_file: fields/simplest_cube_data.msh
  field_name: porosity
  default_value: 0.1
scalar_identic_mesh: !FieldFE
  mesh_data_file: fields/identic_mesh_data.msh
  field_name: scalar
  default_value: 0.0
  interpolation: identic_mesh
vector_identic_mesh: !FieldFE
  mesh_data_file: fields/identic_mesh_data.msh
  field_name: vector_fixed
  default_value: 0.0
  interpolation: identic_mesh
##### tests of intersection interpolation P0
interp_scalar_intersect: !FieldFE
  mesh_data_file: fields/interpolate_boundary_data.msh
  field_name: scalar
  default_value: 0.0
  interpolation: P0_intersection
interp_vector_fixed_intersect: !FieldFE
  mesh_data_file: fields/interpolate_boundary_data.msh
  field_name: vector_fixed
  default_value: 0.0
  interpolation: P0_intersection
##### tests of gauss interpolation P0
interp_scalar_gauss: !FieldFE
  mesh_data_file: fields/interpolate_boundary_data.msh
  field_name: scalar
  default_value: 0.0
interp_vector_fixed_gauss: !FieldFE
  mesh_data_file: fields/interpolate_boundary_data.msh
  field_name: vector_fixed
  default_value: 0.0
#interp_scalar_gauss_large: !FieldFE
#  mesh_data_file: fields/bigger_3d_cube_0.5.msh
#  field_name: scalar
#interp_tensor_gauss_fixed: !FieldFE
#  mesh_data_file: fields/interpolate_boundary_data.msh
#  field_name: tensor_fixed
)YAML";


class FieldFENewTest : public testing::Test {
public:
    typedef FieldFE<3, FieldValue<3>::Scalar > ScalarField;
    typedef FieldFE<3, FieldValue<3>::Enum > EnumField;
    typedef FieldFE<3, FieldValue<3>::VectorFixed > VecFixField;
    typedef FieldFE<3, FieldValue<3>::TensorFixed > TensorField;

    virtual void SetUp() {
        // setup FilePath directories
        FilePath::set_io_dirs(".",UNIT_TESTS_SRC_DIR,"",".");

        Profiler::instance();

        mesh = mesh_full_constructor("{mesh_file=\"mesh/simplest_cube.msh\"}");

        Input::Type::Record rec_type = Input::Type::Record("Test","")
            .declare_key("scalar", ScalarField::get_input_type(), Input::Type::Default::obligatory(),"" )
            .declare_key("scalar_unit_conversion", ScalarField::get_input_type(), Input::Type::Default::obligatory(),"" )
            .declare_key("scalar_time_shift", ScalarField::get_input_type(), Input::Type::Default::obligatory(),"" )
            .declare_key("enum", EnumField::get_input_type(), Input::Type::Default::obligatory(),"" )
            .declare_key("vector_fixed", VecFixField::get_input_type(), Input::Type::Default::obligatory(),"" )
            .declare_key("tensor_fixed", TensorField::get_input_type(), Input::Type::Default::obligatory(),"" )
            .declare_key("vtk_scalar", ScalarField::get_input_type(), Input::Type::Default::obligatory(),"" )
            .declare_key("vtk_vector", VecFixField::get_input_type(), Input::Type::Default::obligatory(),"" )
            .declare_key("vtk_tensor", TensorField::get_input_type(), Input::Type::Default::obligatory(),"" )
            .declare_key("default_values", VecFixField::get_input_type(), Input::Type::Default::obligatory(),"" )
            .declare_key("scalar_identic_mesh", ScalarField::get_input_type(), Input::Type::Default::obligatory(),"" )
            .declare_key("vector_identic_mesh", ScalarField::get_input_type(), Input::Type::Default::obligatory(),"" )
            .declare_key("interp_scalar_intersect", ScalarField::get_input_type(), Input::Type::Default::obligatory(),"" )
            .declare_key("interp_vector_fixed_intersect", VecFixField::get_input_type(), Input::Type::Default::obligatory(),"" )
            .declare_key("interp_scalar_gauss", ScalarField::get_input_type(), Input::Type::Default::obligatory(),"" )
            .declare_key("interp_vector_fixed_gauss", VecFixField::get_input_type(), Input::Type::Default::obligatory(),"" )
//            .declare_key("interp_scalar_gauss_large", ScalarField::get_input_type(), Input::Type::Default::obligatory(),"" )
//            .declare_key("interp_tensor_gauss_fixed", TensorField::get_input_type(), Input::Type::Default::obligatory(),"" )
            .close();

        Input::ReaderToStorage reader( elem_input, rec_type, Input::FileFormat::format_YAML );
        rec=reader.get_root_interface<Input::Record>();

        test_time[0] = 0.0;
        test_time[1] = 1.0;
        test_time[2] = 2.0;

    }
    virtual void TearDown() {
    	delete mesh;
    }

    const FieldAlgoBaseInitData& init_data(std::string field_name) {
    	static const FieldAlgoBaseInitData init_data(field_name, 0, UnitSI::dimensionless());
    	return init_data;
    }

    Mesh * mesh;
    Input::Record rec;
    Space<3>::Point point;
    double test_time[3];

};


TEST_F(FieldFENewTest, scalar) {
    ScalarField field;
    field.init_from_input(rec.val<Input::Record>("scalar"), init_data("scalar"));
    field.set_mesh(mesh,false);

    for (unsigned int j=0; j<2; j++) {
        field.set_time(test_time[j]);
        for(unsigned int i=0; i < mesh->n_elements(); i++) {
            EXPECT_DOUBLE_EQ( j*0.1+(i+1)*0.1 , field.value(point,mesh->element_accessor(i)) );
        }
    }
}


TEST_F(FieldFENewTest, bc_scalar) {
    ScalarField field;
    field.init_from_input(rec.val<Input::Record>("scalar"), init_data("scalar"));
    field.set_mesh(mesh,true);
    for (unsigned int j=0; j<2; j++) {
    	field.set_time(test_time[j]);

        for(unsigned int i=9; i < 15; i++) {
            EXPECT_DOUBLE_EQ( 1.0+j*0.1+(i-8)*0.1 , field.value(point,mesh->element_accessor(i)) );
        }
    }

}


TEST_F(FieldFENewTest, scalar_unit_conv) {
    ScalarField field;
    field.init_from_input(rec.val<Input::Record>("scalar_unit_conversion"), init_data("scalar_unit_conversion"));
    field.set_mesh(mesh,false);

    for (unsigned int j=0; j<2; j++) {
        field.set_time(test_time[j]);
        for(unsigned int i=0; i < mesh->n_elements(); i++) {
            EXPECT_DOUBLE_EQ( j*10.0+(i+1)*10.0 , field.value(point,mesh->element_accessor(i)) );
        }
    }
    field.set_time(test_time[2]); //temporary solution, remove this line after fix 'bc_scalar_unit_conv' test
}


TEST_F(FieldFENewTest, bc_scalar_unit_conv) {
    ScalarField field;
    field.init_from_input(rec.val<Input::Record>("scalar_unit_conversion"), init_data("scalar_unit_conversion"));
    field.set_mesh(mesh,true);
    for (unsigned int j=0; j<3; j++) {
    	field.set_time(test_time[j]);
        for(unsigned int i=9; i < 13; i++) {
            EXPECT_DOUBLE_EQ( 110.0+j*10.0+(i-9)*10.0 , field.value(point,mesh->element_accessor(i)) );
        }
    }

}


TEST_F(FieldFENewTest, scalar_time_shift) {
    ScalarField field;
    field.init_from_input(rec.val<Input::Record>("scalar_time_shift"), init_data("scalar_time_shift"));
    field.set_mesh(mesh,false);

    for (unsigned int j=0; j<2; j++) {
        field.set_time(test_time[j]);
        for(unsigned int i=0; i < mesh->n_elements(); i++) {
            EXPECT_DOUBLE_EQ( j*0.1+(i+2)*0.1 , field.value(point,mesh->element_accessor(i)) );
        }
    }
}


TEST_F(FieldFENewTest, vector_fixed) {
	string expected_vals[2] = {"1 2 3", "2 3 4"};
    VecFixField field;
    field.init_from_input(rec.val<Input::Record>("vector_fixed"), init_data("vector_fixed"));
    field.set_mesh(mesh,false);
     for (unsigned int j=0; j<2; j++) {
    	field.set_time(test_time[j]);
         for(unsigned int i=0; i < mesh->n_elements(); i++) {
            EXPECT_TRUE( arma::min(arma::vec3(expected_vals[j]) == field.value(point,mesh->element_accessor(i))) );
        }
    }
}


TEST_F(FieldFENewTest, bc_vector_fixed) {
	string expected_vals[2] = {"4 5 6", "5 6 7"};
    VecFixField field;
    field.init_from_input(rec.val<Input::Record>("vector_fixed"), init_data("vector_fixed"));
    field.set_mesh(mesh,true);
     for (unsigned int j=0; j<2; j++) {
    	field.set_time(test_time[j]);
     	for(unsigned int i=9; i < 15; i++) {
            EXPECT_TRUE( arma::min(arma::vec3(expected_vals[j]) == field.value(point,mesh->element_accessor(i))) );
        }
    }
}


TEST_F(FieldFENewTest, tensor_fixed) {
	string expected_vals[2] = {"1 2 3; 4 5 6; 7 8 9", "2 3 4; 5 6 7; 8 9 10"};
    TensorField field;
    field.init_from_input(rec.val<Input::Record>("tensor_fixed"), init_data("tensor_fixed"));
    field.set_mesh(mesh,false);
     for (unsigned int j=0; j<2; j++) {
    	field.set_time(test_time[j]);
     	for(unsigned int i=0; i < mesh->n_elements(); i++) {
    		arma::umat match = ( arma::mat33(expected_vals[j]) == field.value(point,mesh->element_accessor(i)) );
            EXPECT_TRUE( match.min() );
        }
    }
}


TEST_F(FieldFENewTest, bc_tensor_fixed) {
	string expected_vals[2] = {"4 5 6; 7 8 9; 10 11 12", "5 6 7; 8 9 10; 11 12 13"};
    TensorField field;
    field.init_from_input(rec.val<Input::Record>("tensor_fixed"), init_data("tensor_fixed"));
    field.set_mesh(mesh, true);
     for (unsigned int j=0; j<2; j++) {
    	field.set_time(test_time[j]);
         for(unsigned int i=9; i < 15; i++) {
            arma::umat match = ( arma::mat33(expected_vals[j]) == field.value(point,mesh->element_accessor(i)) );
            EXPECT_TRUE( match.min() );
        }
    }
}


TEST_F(FieldFENewTest, vtk_scalar) {
	ScalarField field;
    field.init_from_input(rec.val<Input::Record>("vtk_scalar"), init_data("vtk_scalar"));
    field.set_mesh(mesh, false);
	field.set_time(0.0);

	for(unsigned int i=0; i<mesh->n_elements(); i++) {
		EXPECT_DOUBLE_EQ( (i+1)*0.1 , field.value(point,mesh->element_accessor(i)) );
	}
}



TEST_F(FieldFENewTest, vtk_vector) {
	string expected_vals = "0.5 1 1.5";
    VecFixField field;
    field.init_from_input(rec.val<Input::Record>("vtk_vector"), init_data("vtk_vector"));
    field.set_mesh(mesh, false);
   	field.set_time(0.0);
    for(unsigned int i=0; i < mesh->n_elements(); i++) {
    	EXPECT_TRUE( arma::min(arma::vec3(expected_vals) == field.value(point,mesh->element_accessor(i))) );
    }
}


TEST_F(FieldFENewTest, vtk_tensor) {
	string expected_vals = "1 2 3; 4 5 6; 7 8 9";
    TensorField field;
    field.init_from_input(rec.val<Input::Record>("vtk_tensor"), init_data("vtk_tensor"));
    field.set_mesh(mesh,false);
   	field.set_time(0.0);
    for(unsigned int i=0; i < mesh->n_elements(); i++) {
    	arma::umat match = ( arma::mat33(expected_vals) == field.value(point,mesh->element_accessor(i)) );
        EXPECT_TRUE( match.min() );
    }
}


TEST_F(FieldFENewTest, scalar_enum) {
    EnumField field;
    field.init_from_input(rec.val<Input::Record>("enum"), init_data("enum"));
    field.set_mesh(mesh,false);
    for (unsigned int j=0; j<2; j++) {
    	field.set_time(test_time[j]);
     	for(unsigned int i=0; i < mesh->n_elements(); i++) {
            EXPECT_EQ( j, field.value(point,mesh->element_accessor(i)) );
        }
    }
}


TEST_F(FieldFENewTest, bc_scalar_enum) {
    EnumField field;
    field.init_from_input(rec.val<Input::Record>("enum"), init_data("enum"));
    field.set_mesh(mesh, true);
    for (unsigned int j=0; j<2; j++) {
		field.set_time(test_time[j]);
 		for(unsigned int i=9; i < 13; i++) {
			EXPECT_EQ( j+1, field.value(point,mesh->element_accessor(i)) );
		}
    }
}


TEST_F(FieldFENewTest, default_values) {
	string expected_vals = "0.1 0.1 0.1";
    VecFixField field;
    field.init_from_input(rec.val<Input::Record>("default_values"), init_data("default_values"));
    field.set_mesh(mesh,true);
    for (unsigned int j=0; j<2; j++) {
    	field.set_time(test_time[j]);
     	for(unsigned int i=9; i < 13; i++) {
            EXPECT_TRUE( arma::min(arma::vec3(expected_vals) == field.value(point,mesh->element_accessor(i))) );
        }
    }
}


TEST_F(FieldFENewTest, scalar_identic_mesh) {
    ScalarField field;
    field.init_from_input(rec.val<Input::Record>("scalar_identic_mesh"), init_data("scalar_identic_mesh"));
    field.set_mesh(mesh,false);

    for (unsigned int j=0; j<2; j++) {
        field.set_time(test_time[j]);
        for(unsigned int i=0; i < mesh->n_elements(); i++) {
            EXPECT_DOUBLE_EQ( 1.0+j*0.1+(i+1)*0.1 , field.value(point,mesh->element_accessor(i)) );
        }
    }
}


TEST_F(FieldFENewTest, bc_scalar_identic_mesh) {
    ScalarField field;
    field.init_from_input(rec.val<Input::Record>("scalar_identic_mesh"), init_data("scalar_identic_mesh"));
    field.set_mesh(mesh,true);
    for (unsigned int j=0; j<2; j++) {
    	field.set_time(test_time[j]);

        for(unsigned int i=9; i < 13; i++) {
            EXPECT_DOUBLE_EQ( 2.0+j*0.1+(i-8)*0.1 , field.value(point,mesh->element_accessor(i)) );
        }
    }

}


TEST_F(FieldFENewTest, vector_fixed_identic_mesh) {
	string expected_vals[2] = {"3 4 5", "6 7 8"};
    VecFixField field;
    field.init_from_input(rec.val<Input::Record>("vector_identic_mesh"), init_data("vector_identic_mesh"));
    field.set_mesh(mesh,false);
     for (unsigned int j=0; j<2; j++) {
    	field.set_time(test_time[j]);
         for(unsigned int i=0; i < mesh->n_elements(); i++) {
            EXPECT_TRUE( arma::min(arma::vec3(expected_vals[j]) == field.value(point,mesh->element_accessor(i))) );
        }
    }
}


TEST_F(FieldFENewTest, bc_vector_fixed_identic_mesh) {
	string expected_vals[2] = {"1 2 3", "4 5 6"};
    VecFixField field;
    field.init_from_input(rec.val<Input::Record>("vector_identic_mesh"), init_data("vector_identic_mesh"));
    field.set_mesh(mesh,true);
     for (unsigned int j=0; j<2; j++) {
    	field.set_time(test_time[j]);
     	for(unsigned int i=9; i < 13; i++) {
            EXPECT_TRUE( arma::min(arma::vec3(expected_vals[j]) == field.value(point,mesh->element_accessor(i))) );
        }
    }
}


TEST_F(FieldFENewTest, intersection_1d_2d_elements_small_scalar) {
    ScalarField field;
    field.init_from_input(rec.val<Input::Record>("interp_scalar_intersect"), init_data("interp_scalar_intersect"));
    field.set_mesh(mesh, true);
    //std::vector<unsigned int> expected_vals = {4,9,6,7};

    for (unsigned int j=0; j<2; j++) {
    	field.set_time(test_time[j]);
    	std::cout << "Time: " << test_time[j] << std::endl;

    	for (unsigned int i=9; i<13; ++i) {
    		ElementAccessor<3> elm = mesh->element_accessor(i);
    		std::cout << " - " << field.value(elm.centre(), elm) << std::endl;
    		//EXPECT_DOUBLE_EQ( 0.1*(j+expected_vals[i]), field.value(point, mesh->element_accessor(i+9)) );
    	}
    }

} // */


TEST_F(FieldFENewTest, intersection_1d_2d_elements_small_vector) {
	VecFixField field;
    field.init_from_input(rec.val<Input::Record>("interp_vector_fixed_intersect"), init_data("interp_vector_fixed_intersect"));
    field.set_mesh(mesh, true);
    //std::vector<unsigned int> expected_vals = {4,9,6,7};

    for (unsigned int j=0; j<2; j++) {
    	field.set_time(test_time[j]);
    	std::cout << "Time: " << test_time[j] << std::endl;

    	for (unsigned int i=9; i<13; ++i) {
    		ElementAccessor<3> elm = mesh->element_accessor(i);
    		std::cout << " - " << field.value(elm.centre(), elm) << std::endl;
    		//EXPECT_DOUBLE_EQ( 0.1*(j+expected_vals[i]), field.value(point, mesh->element_accessor(i+9)) );
    	}
    }

} // */


TEST_F(FieldFENewTest, gauss_1d_2d_elements_small_scalar) {
    ScalarField field;
    field.init_from_input(rec.val<Input::Record>("interp_scalar_gauss"), init_data("interp_scalar_gauss"));
    field.set_mesh(mesh, true);
    std::vector<unsigned int> expected_vals = {4,9,6,7};

    for (unsigned int j=0; j<2; j++) {
    	field.set_time(test_time[j]);
    	std::cout << "Time: " << test_time[j] << std::endl;

    	for (unsigned int i=0; i<4; ++i) {
    		ElementAccessor<3> elm = mesh->element_accessor(i+9);
    		std::cout << " - " << field.value(elm.centre(), elm) << std::endl;
    		//EXPECT_DOUBLE_EQ( 0.1*(j+expected_vals[i]), field.value(point, mesh->element_accessor(i+9)) );
    	}
    }

} // */


TEST_F(FieldFENewTest, gauss_1d_2d_elements_small_vector) {
	VecFixField field;
    field.init_from_input(rec.val<Input::Record>("interp_vector_fixed_gauss"), init_data("interp_vector_fixed_gauss"));
    field.set_mesh(mesh, true);
    std::vector<unsigned int> expected_vals = {4,9,6,7};

    for (unsigned int j=0; j<2; j++) {
    	field.set_time(test_time[j]);
    	std::cout << "Time: " << test_time[j] << std::endl;

    	for (unsigned int i=0; i<4; ++i) {
    		ElementAccessor<3> elm = mesh->element_accessor(i+9);
    		std::cout << " - " << field.value(elm.centre(), elm) << std::endl;
    		//EXPECT_DOUBLE_EQ( 0.1*(j+expected_vals[i]), field.value(point, mesh->element_accessor(i+9)) );
    	}
    }

} // */
