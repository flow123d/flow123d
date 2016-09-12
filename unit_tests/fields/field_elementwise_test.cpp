/*
 * field_elementwise_test.cpp
 *
 *  Created on: Jan 25, 2013
 *      Author: jb
 */



#define TEST_USE_PETSC
#define FEAL_OVERRIDE_ASSERTS
#include <flow_gtest_mpi.hh>


#include "fields/field_elementwise.hh"
#include "input/input_type.hh"
#include "input/accessors.hh"
#include "input/reader_to_storage.hh"

#include "system/sys_profiler.hh"

#include "mesh/mesh.h"
#include "mesh/msh_gmshreader.h"


string input = R"INPUT(
{   
   scalar={
       TYPE="FieldElementwise",
       gmsh_file="fields/simplest_cube_data.msh",
       field_name="scalar"
   },
   vector_fixed={
       TYPE="FieldElementwise",
       gmsh_file="fields/simplest_cube_data.msh",
       field_name="vector_fixed"
   },
   tensor_fixed={
       TYPE="FieldElementwise",
       gmsh_file="fields/simplest_cube_data.msh",
       field_name="tensor_fixed"
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

        Profiler::initialize();
        
        FilePath mesh_file( "mesh/simplest_cube.msh", FilePath::input_file);
        mesh= new Mesh;
        ifstream in(string( mesh_file ).c_str());
        mesh->read_gmsh_from_stream(in);

        Input::Type::Record rec_type = Input::Type::Record("Test","")
            .declare_key("scalar", ScalarField::get_input_type(), Input::Type::Default::obligatory(),"" )
            .declare_key("vector_fixed", VecFixField::get_input_type(), Input::Type::Default::obligatory(),"" )
            .declare_key("tensor_fixed", TensorField::get_input_type(), Input::Type::Default::obligatory(),"" )
            .close();

        Input::ReaderToStorage reader( input, rec_type, Input::FileFormat::format_JSON );
        rec=reader.get_root_interface<Input::Record>();

        test_time[0] = 0.0;
        test_time[1] = 1.0;

    }
    virtual void TearDown() {

    }

    Mesh *mesh;
    Input::Record rec;
    Space<3>::Point point;
    double test_time[2];

};


TEST_F(FieldElementwiseTest, scalar) {
    ScalarField field;
    field.init_from_input(rec.val<Input::Record>("scalar"));
    field.set_mesh(mesh,false);

    for (unsigned int j=0; j<2; j++) {
        field.set_time(test_time[j]);
        for(unsigned int i=0; i < mesh->element.size(); i++) {
            EXPECT_DOUBLE_EQ( j*0.1+(i+1)*0.1 , field.value(point,mesh->element_accessor(i)) );
        }
    }
}



TEST_F(FieldElementwiseTest, bc_scalar) {
    ScalarField field;
    field.set_mesh(mesh,true);
    field.init_from_input(rec.val<Input::Record>("scalar"));

    for (unsigned int j=0; j<2; j++) {
    	field.set_time(test_time[j]);

        for(unsigned int i=0; i < 4; i++) {
            EXPECT_DOUBLE_EQ( 1.0+j*0.1+(i+1)*0.1 , field.value(point,mesh->element_accessor(i, true)) );
        }
        EXPECT_DOUBLE_EQ( 0.0, field.value(point,mesh->element_accessor(5, true)) );
    }

}

TEST_F(FieldElementwiseTest, vector_fixed) {
	string expected_vals[2] = {"1 2 3", "2 3 4"};
    VecFixField field;
    field.init_from_input(rec.val<Input::Record>("vector_fixed"));
    field.set_mesh(mesh,false);

    for (unsigned int j=0; j<2; j++) {
    	field.set_time(test_time[j]);

        for(unsigned int i=0; i < mesh->element.size(); i++) {
            EXPECT_TRUE( arma::min(arma::vec3(expected_vals[j]) == field.value(point,mesh->element_accessor(i))) );
        }
    }
}



TEST_F(FieldElementwiseTest, bc_vector_fixed) {
	string expected_vals[2] = {"4 5 6", "5 6 7"};
    VecFixField field;
    field.init_from_input(rec.val<Input::Record>("vector_fixed"));
    field.set_mesh(mesh,true);

    for (unsigned int j=0; j<2; j++) {
    	field.set_time(test_time[j]);

    	for(unsigned int i=0; i < 4; i++) {
            EXPECT_TRUE( arma::min(arma::vec3(expected_vals[j]) == field.value(point,mesh->element_accessor(i,true))) );
        }
        EXPECT_TRUE( arma::min(arma::vec3("0 0 0") == field.value(point,mesh->element_accessor(5,true))) );
    }
}

TEST_F(FieldElementwiseTest, tensor_fixed) {
	string expected_vals[2] = {"1 4 7; 2 5 8; 3 6 9", "2 5 8; 3 6 9; 4 7 10"};
    TensorField field;
    field.init_from_input(rec.val<Input::Record>("tensor_fixed"));
    field.set_mesh(mesh,false);

    for (unsigned int j=0; j<2; j++) {
    	field.set_time(test_time[j]);

    	for(unsigned int i=0; i < mesh->element.size(); i++) {
    		arma::umat match = ( arma::mat33(expected_vals[j]) == field.value(point,mesh->element_accessor(i)) );
            EXPECT_TRUE( match.min() );
        }
    }
}




TEST_F(FieldElementwiseTest, bc_tensor_fixed) {
	string expected_vals[2] = {"4 7 10; 5 8 11; 6 9 12", "5 8 11; 6 9 12; 7 10 13"};
    TensorField field;
    field.init_from_input(rec.val<Input::Record>("tensor_fixed"));
    field.set_mesh(mesh, true);

    for (unsigned int j=0; j<2; j++) {
    	field.set_time(test_time[j]);

        for(unsigned int i=0; i < 4; i++) {
            arma::umat match = ( arma::mat33(expected_vals[j]) == field.value(point,mesh->element_accessor(i,true)) );
            EXPECT_TRUE( match.min() );
        }
        arma::umat match = ( arma::mat33("0 0 0; 0 0 0; 0 0 0") == field.value(point,mesh->element_accessor(5,true)) );
        EXPECT_TRUE( match.min() );
    }
}



TEST_F(FieldElementwiseTest, scalar_enum) {
    EnumField field;
    field.set_mesh(mesh,false);

    for (unsigned int j=0; j<2; j++) {
    	field.set_time(test_time[j]);

    	for(unsigned int i=0; i < mesh->element.size(); i++) {
            EXPECT_EQ( (unsigned int)0, field.value(point,mesh->element_accessor(i)) );
        }
    }
}




TEST_F(FieldElementwiseTest, bc_scalar_enum) {
    EnumField field;
    field.set_mesh(mesh, true);

    for (unsigned int j=0; j<2; j++) {
		field.set_time(test_time[j]);

		for(unsigned int i=0; i<6; i++) {
			unsigned int val = i + j + ( i<4 ? 1 : 10 );
			field.set_data_row(i, val );
		}
		for(unsigned int i=0; i < 4; i++) {
			EXPECT_EQ( i+j+1, field.value(point,mesh->element_accessor(i,true)) );
		}
		EXPECT_EQ( 14+j, field.value(point,mesh->element_accessor(4,true)) );
		EXPECT_EQ( 15+j, field.value(point,mesh->element_accessor(5,true)) );
    }
}

