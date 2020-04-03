/*
 * field_elementwise_test.cpp
 *
 *  Created on: Jan 25, 2013
 *      Author: jb
 */



#define TEST_USE_PETSC
#define FEAL_OVERRIDE_ASSERTS
#include <flow_gtest_mpi.hh>
#include <mesh_constructor.hh>


#include "fields/field_elementwise.hh"
#include "tools/unit_si.hh"
#include "input/input_type.hh"
#include "input/accessors.hh"
#include "input/reader_to_storage.hh"

#include "system/sys_profiler.hh"

#include "mesh/mesh.h"
#include "io/msh_gmshreader.h"



string input = R"INPUT(
{   
   scalar={
       TYPE="FieldElementwise",
       mesh_data_file="fields/simplest_cube_data.msh",
       field_name="scalar",
       default_value=0.0
   },
   scalar_unit_conversion={
       TYPE="FieldElementwise",
       mesh_data_file="fields/simplest_cube_data.msh",
       field_name="scalar",
       unit="const; const=100*m^0",
       default_value=0.0
   },
   scalar_time_shift={
       TYPE="FieldElementwise",
       mesh_data_file="fields/simplest_cube_data.msh",
       field_name="scalar",
       default_value=0.0,
       read_time_shift=1.0
   },
   vector_fixed={
       TYPE="FieldElementwise",
       mesh_data_file="fields/simplest_cube_data.msh",
       field_name="vector_fixed",
       default_value=0.0
   },
   tensor_fixed={
       TYPE="FieldElementwise",
       mesh_data_file="fields/simplest_cube_data.msh",
       field_name="tensor_fixed",
       default_value=0.0
   },
   vtk_scalar={
       TYPE="FieldElementwise",
       mesh_data_file="fields/vtk_ascii_data.vtu",
       field_name="scalar_field"
   },
   vtk_vector={
       TYPE="FieldElementwise",
       mesh_data_file="fields/vtk_binary_data.vtu",
       field_name="vector_field"
   },
   vtk_tensor={
       TYPE="FieldElementwise",
       mesh_data_file="fields/vtk_compressed_data.vtu",
       field_name="tensor_field"
   },
   default_values={
       TYPE="FieldElementwise",
       mesh_data_file="fields/simplest_cube_data.msh",
       field_name="porosity",
       default_value=0.1
   }
}
)INPUT";


class FieldElementwiseTest : public testing::Test {
public:
    typedef FieldElementwise<3, FieldValue<3>::Scalar > ScalarField;
    typedef FieldElementwise<3, FieldValue<3>::Enum > EnumField;
    typedef FieldElementwise<3, FieldValue<3>::VectorFixed > VecFixField;
    typedef FieldElementwise<3, FieldValue<3>::TensorFixed > TensorField;

    virtual void SetUp() {
        // setup FilePath directories
        FilePath::set_io_dirs(".",UNIT_TESTS_SRC_DIR,"",".");

        Profiler::instance();
        
        mesh = mesh_full_constructor("{mesh_file=\"mesh/simplest_cube.msh\"}");

        Input::Type::Record rec_type = Input::Type::Record("Test","")
            .declare_key("scalar", ScalarField::get_input_type(), Input::Type::Default::obligatory(),"" )
            .declare_key("scalar_unit_conversion", ScalarField::get_input_type(), Input::Type::Default::obligatory(),"" )
            .declare_key("scalar_time_shift", ScalarField::get_input_type(), Input::Type::Default::obligatory(),"" )
            .declare_key("vector_fixed", VecFixField::get_input_type(), Input::Type::Default::obligatory(),"" )
            .declare_key("tensor_fixed", TensorField::get_input_type(), Input::Type::Default::obligatory(),"" )
            .declare_key("vtk_scalar", ScalarField::get_input_type(), Input::Type::Default::obligatory(),"" )
            .declare_key("vtk_vector", VecFixField::get_input_type(), Input::Type::Default::obligatory(),"" )
            .declare_key("vtk_tensor", TensorField::get_input_type(), Input::Type::Default::obligatory(),"" )
            .declare_key("default_values", VecFixField::get_input_type(), Input::Type::Default::obligatory(),"" )
            .close();

        Input::ReaderToStorage reader( input, rec_type, Input::FileFormat::format_JSON );
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


TEST_F(FieldElementwiseTest, scalar) {
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



TEST_F(FieldElementwiseTest, bc_scalar) {
    ScalarField field;
    field.init_from_input(rec.val<Input::Record>("scalar"), init_data("scalar"));
    field.set_mesh(mesh,true);

    for (unsigned int j=0; j<2; j++) {
    	field.set_time(test_time[j]);

        for(unsigned int i=9; i < 13; i++) {
            EXPECT_DOUBLE_EQ( 1.0+j*0.1+(i-8)*0.1 , field.value(point,mesh->element_accessor(i)) );
        }
        EXPECT_DOUBLE_EQ( 0.0, field.value(point,mesh->element_accessor(14)) );
    }

}

TEST_F(FieldElementwiseTest, scalar_unit_conv) {
    ScalarField field;
    field.init_from_input(rec.val<Input::Record>("scalar_unit_conversion"), init_data("scalar_unit_conversion"));
    field.set_mesh(mesh,false);

    for (unsigned int j=0; j<2; j++) {
        field.set_time(test_time[j]);
        for(unsigned int i=0; i < mesh->n_elements(); i++) {
            EXPECT_DOUBLE_EQ( j*10.0+(i+1)*10.0 , field.value(point,mesh->element_accessor(i)) );
        }
    }
}

TEST_F(FieldElementwiseTest, bc_scalar_unit_conv) {
    ScalarField field;
    field.init_from_input(rec.val<Input::Record>("scalar_unit_conversion"), init_data("scalar_unit_conversion"));
    field.set_mesh(mesh,true);

    for (unsigned int j=0; j<3; j++) {
    	field.set_time(test_time[j]);

        for(unsigned int i=9; i < 13; i++) {
            EXPECT_DOUBLE_EQ( 110.0+j*10.0+(i-9)*10.0 , field.value(point,mesh->element_accessor(i)) );
        }
        EXPECT_DOUBLE_EQ( 0.0, field.value(point,mesh->element_accessor(13)) );
    }

}

TEST_F(FieldElementwiseTest, scalar_time_shift) {
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

TEST_F(FieldElementwiseTest, vector_fixed) {
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



TEST_F(FieldElementwiseTest, bc_vector_fixed) {
	string expected_vals[2] = {"4 5 6", "5 6 7"};
    VecFixField field;
    field.init_from_input(rec.val<Input::Record>("vector_fixed"), init_data("vector_fixed"));
    field.set_mesh(mesh,true);

    for (unsigned int j=0; j<2; j++) {
    	field.set_time(test_time[j]);

    	for(unsigned int i=9; i < 13; i++) {
            EXPECT_TRUE( arma::min(arma::vec3(expected_vals[j]) == field.value(point,mesh->element_accessor(i))) );
        }
        EXPECT_TRUE( arma::min(arma::vec3("0 0 0") == field.value(point,mesh->element_accessor(13))) );
    }
}

TEST_F(FieldElementwiseTest, tensor_fixed) {
	string expected_vals[2] = {"1 4 7; 2 5 8; 3 6 9", "2 5 8; 3 6 9; 4 7 10"};
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




TEST_F(FieldElementwiseTest, bc_tensor_fixed) {
	string expected_vals[2] = {"4 7 10; 5 8 11; 6 9 12", "5 8 11; 6 9 12; 7 10 13"};
    TensorField field;
    field.init_from_input(rec.val<Input::Record>("tensor_fixed"), init_data("tensor_fixed"));
    field.set_mesh(mesh, true);

    for (unsigned int j=0; j<2; j++) {
    	field.set_time(test_time[j]);

        for(unsigned int i=9; i < 13; i++) {
            arma::umat match = ( arma::mat33(expected_vals[j]) == field.value(point,mesh->element_accessor(i)) );
            EXPECT_TRUE( match.min() );
        }
        arma::umat match = ( arma::mat33("0 0 0; 0 0 0; 0 0 0") == field.value(point,mesh->element_accessor(13)) );
        EXPECT_TRUE( match.min() );
    }
}



TEST_F(FieldElementwiseTest, vtk_scalar) {
	ScalarField field;
    field.init_from_input(rec.val<Input::Record>("vtk_scalar"), init_data("vtk_scalar"));
    field.set_mesh(mesh, false);
	field.set_time(0.0);

	for(unsigned int i=0; i<mesh->n_elements(); i++) {
		EXPECT_DOUBLE_EQ( (i+1)*0.1 , field.value(point,mesh->element_accessor(i)) );
	}
}



TEST_F(FieldElementwiseTest, vtk_vector) {
	string expected_vals = "0.5 1 1.5";
    VecFixField field;
    field.init_from_input(rec.val<Input::Record>("vtk_vector"), init_data("vtk_vector"));
    field.set_mesh(mesh, false);
   	field.set_time(0.0);

    for(unsigned int i=0; i < mesh->n_elements(); i++) {
    	EXPECT_TRUE( arma::min(arma::vec3(expected_vals) == field.value(point,mesh->element_accessor(i))) );
    }
}



TEST_F(FieldElementwiseTest, vtk_tensor) {
	string expected_vals = "1 4 7; 2 5 8; 3 6 9";
    TensorField field;
    field.init_from_input(rec.val<Input::Record>("vtk_tensor"), init_data("vtk_tensor"));
    field.set_mesh(mesh,false);
   	field.set_time(0.0);

    for(unsigned int i=0; i < mesh->n_elements(); i++) {
    	arma::umat match = ( arma::mat33(expected_vals) == field.value(point,mesh->element_accessor(i)) );
        EXPECT_TRUE( match.min() );
    }
}


TEST_F(FieldElementwiseTest, scalar_enum) {
    EnumField field;
    field.set_mesh(mesh,false);

    for (unsigned int j=0; j<2; j++) {
    	field.set_time(test_time[j]);

    	for(unsigned int i=0; i < mesh->n_elements(); i++) {
            EXPECT_EQ( (unsigned int)0, field.value(point,mesh->element_accessor(i)) );
        }
    }
}




TEST_F(FieldElementwiseTest, bc_scalar_enum) {
    EnumField field;
    field.set_mesh(mesh, true);

    for (unsigned int j=0; j<2; j++) {
		field.set_time(test_time[j]);

		//for(unsigned int i=0; i<6; i++) {
		//	unsigned int val = i + j + ( i<4 ? 1 : 10 );
		//	field.set_data_row(i, val );
		//}
		for(unsigned int i=9; i < 15; i++) {
			EXPECT_EQ( (unsigned int)0, field.value(point,mesh->element_accessor(i)) );
		}
		//EXPECT_EQ( 14+j, field.value(point,mesh->element_accessor(4)) );
		//EXPECT_EQ( 15+j, field.value(point,mesh->element_accessor(5)) );
    }
}




TEST_F(FieldElementwiseTest, default_values) {
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
