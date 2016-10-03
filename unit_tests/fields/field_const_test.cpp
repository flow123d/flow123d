/*
 * field_const_test.cpp
 *
 *  Created on: Dec 18, 2012
 *      Author: jb
 */



#define FEAL_OVERRIDE_ASSERTS

#include <flow_gtest.hh>
#include <memory>


#include "fields/field_constant.hh"
#include "fields/table_function.hh"
#include "input/input_type.hh"
#include "input/accessors.hh"
#include "input/reader_to_storage.hh"


FLOW123D_FORCE_LINK_IN_PARENT(field_constant)


string field_const_input = R"YAML(
tensor1: 3.14
tensor2:
 - 1
 - 2
 - 3
tensor3:
 - 1
 - 2
 - 3
 - 4
 - 5
 - 6
tensor4: [ [1,2,3], [4,5,6], [7,8,9]]
init_conc: !FieldConstant
 - value:
   - 1.2
   - 2.3
   - 3.4
)YAML";


typedef FieldAlgorithmBase<3, FieldValue<3>::TensorFixed > TensorField;
typedef FieldAlgorithmBase<3, FieldValue<3>::VectorFixed > VectorField;
typedef FieldAlgorithmBase<3, FieldValue<0>::Scalar > ScalarField;


void check_tensor_field(std::shared_ptr<TensorField> field,const  std::string& expected, vector<Space<3>::Point> points,ElementAccessor<3> elem) {
    arma::mat::fixed<3,3> result;
    arma::mat::fixed<3,3> exp(expected);

    for(auto point : points) {
    	result = field->value( point, elem);
    	arma::umat mat_match = (  exp == result );
    	bool mat_equal = arma::min( arma::min(mat_match));
    	if (! mat_equal) cout << "Expected: " << exp << endl << "Result: "<< result << endl;
    	EXPECT_TRUE( mat_equal );
    }
}


TEST(FieldConst, read_from_input) {

    // setup FilePath directories
    FilePath::set_io_dirs(".",UNIT_TESTS_SRC_DIR,"",".");

    Input::Type::Record rec_type = Input::Type::Record("FieldConstTest","")
        .declare_key("tensor1", TensorField::get_input_type_instance(), Input::Type::Default::obligatory(),"" )
        .declare_key("tensor2", TensorField::get_input_type_instance(), Input::Type::Default::obligatory(),"" )
        .declare_key("tensor3", TensorField::get_input_type_instance(), Input::Type::Default::obligatory(),"" )
        .declare_key("tensor4", TensorField::get_input_type_instance(), Input::Type::Default::obligatory(),"" )
        //.declare_key("init_conc", VectorField::get_input_type_instance(), Input::Type::Default::obligatory(), "" )
        .close();

    // read input string
    Input::ReaderToStorage reader( field_const_input, rec_type, Input::FileFormat::format_YAML );
    Input::Record in_rec=reader.get_root_interface<Input::Record>();

    Space<3>::Point point_1, point_2;
    point_1(0)=1.0; point_1(1)=2.0; point_1(2)=3.0;
    point_2(0)=2.0; point_2(1)=4.0; point_2(2)=6.0;
    ElementAccessor<3> elm;

    /*
    auto conc=VectorField::function_factory(in_rec.val<Input::AbstractRecord>("init_conc"), 3);
    {
        arma::vec result;

        result = conc->value( point_1, elm);
        EXPECT_DOUBLE_EQ( 1.2 , result[0]);
        EXPECT_DOUBLE_EQ( 2.3, result[1]);
        EXPECT_DOUBLE_EQ( 3.4, result[2]);

        result = conc->value( point_2, elm);
        EXPECT_DOUBLE_EQ( 1.2 , result[0]);
        EXPECT_DOUBLE_EQ( 2.3, result[1]);
        EXPECT_DOUBLE_EQ( 3.4, result[2]);
    }*/

    auto tensor1=TensorField::function_factory(in_rec.val<Input::AbstractRecord>("tensor1"));
    check_tensor_field(tensor1, "3.14 0 0; 0 3.14 0; 0 0 3.14", {point_1, point_2}, elm);

    auto tensor2=TensorField::function_factory(in_rec.val<Input::AbstractRecord>("tensor2"));
    check_tensor_field(tensor2, "1 0 0; 0 2 0; 0 0 3", {point_1, point_2}, elm);

    auto tensor3=TensorField::function_factory(in_rec.val<Input::AbstractRecord>("tensor3"));
    check_tensor_field(tensor3, "1 2 3; 2 4 5; 3 5 6", {point_1, point_2}, elm);

    auto tensor4=TensorField::function_factory(in_rec.val<Input::AbstractRecord>("tensor4"));
    check_tensor_field(tensor4, "1 2 3; 4 5 6; 7 8 9", {point_1, point_2}, elm);
}


string field_const_table_function_input = R"YAML(
table_function_scalar: !FieldConstant
  value: 0.5
  time_function:
    - - 0.0
      - 0.5
    - - 1.0
      - 1.0
    - - 2.0
      - 3.0
table_function_vector: !FieldConstant
  value: [0.5, 1.5, 2.0]
  time_function:
    - - 0.0
      - [0.5, 1.5, 2.0]
    - - 1.0
      - [1.0, 1.5, 1.0]
    - - 2.0
      - [3.0, 1.5, 5.0]
table_function_tensor: !FieldConstant
  value: [ [1,3,4], [0,3,4], [1,6,6] ]
  time_function:
    - - 0.0
      - [ [1,3,4], [0,3,4], [1,6,6] ]
    - - 1.0
      - [ [3,3,6], [0,2,5], [2,7,4] ]
    - - 2.0
      - [ [5,3,4], [2,3,7], [5,6,3] ]
)YAML";

TEST(FieldConst, table_function) {

    // setup FilePath directories
    FilePath::set_io_dirs(".",UNIT_TESTS_SRC_DIR,"",".");

    Input::Type::Record rec_type = Input::Type::Record("FieldConstTableFunctionTest","")
        .declare_key("table_function_scalar", ScalarField::get_input_type_instance(), Input::Type::Default::obligatory(),"" )
        .declare_key("table_function_vector", VectorField::get_input_type_instance(), Input::Type::Default::obligatory(),"" )
        .declare_key("table_function_tensor", TensorField::get_input_type_instance(), Input::Type::Default::obligatory(),"" )
        .close();

    // read input string
    Input::ReaderToStorage reader( field_const_table_function_input, rec_type, Input::FileFormat::format_YAML );
    Input::Record in_rec=reader.get_root_interface<Input::Record>();

    {
        std::shared_ptr< ScalarField > scalar_base = ScalarField::function_factory(in_rec.val<Input::AbstractRecord>("table_function_scalar"));
        auto scalar = std::static_pointer_cast< FieldConstant<3, FieldValue<0>::Scalar> >(scalar_base);

        EXPECT_DOUBLE_EQ( scalar->time_dependent_value(0.0), 0.5 );
    	EXPECT_DOUBLE_EQ( scalar->time_dependent_value(0.2), 0.6 );
    	EXPECT_DOUBLE_EQ( scalar->time_dependent_value(0.4), 0.7 );
    	EXPECT_DOUBLE_EQ( scalar->time_dependent_value(0.6), 0.8 );
    	EXPECT_DOUBLE_EQ( scalar->time_dependent_value(0.8), 0.9 );
    	EXPECT_DOUBLE_EQ( scalar->time_dependent_value(1.0), 1.0 );
    	EXPECT_DOUBLE_EQ( scalar->time_dependent_value(1.5), 2.0 );
    	EXPECT_DOUBLE_EQ( scalar->time_dependent_value(2.0), 3.0 );
    }

    {
        std::shared_ptr< VectorField > vector_base = VectorField::function_factory(in_rec.val<Input::AbstractRecord>("table_function_vector"));
        auto vector = std::static_pointer_cast< FieldConstant<3, FieldValue<3>::VectorFixed> >(vector_base);
    	arma::vec result;

    	result = vector->time_dependent_value(0.2);
        EXPECT_DOUBLE_EQ( 0.6, result[0]);
        EXPECT_DOUBLE_EQ( 1.5, result[1]);
        EXPECT_DOUBLE_EQ( 1.8, result[2]);

    	result = vector->time_dependent_value(0.6);
        EXPECT_DOUBLE_EQ( 0.8, result[0]);
        EXPECT_DOUBLE_EQ( 1.5, result[1]);
        EXPECT_DOUBLE_EQ( 1.4, result[2]);

    	result = vector->time_dependent_value(1.0);
        EXPECT_DOUBLE_EQ( 1.0, result[0]);
        EXPECT_DOUBLE_EQ( 1.5, result[1]);
        EXPECT_DOUBLE_EQ( 1.0, result[2]);

    	result = vector->time_dependent_value(1.5);
        EXPECT_DOUBLE_EQ( 2.0, result[0]);
        EXPECT_DOUBLE_EQ( 1.5, result[1]);
        EXPECT_DOUBLE_EQ( 3.0, result[2]);

    }

    {
        std::shared_ptr< TensorField > tensor_base = TensorField::function_factory(in_rec.val<Input::AbstractRecord>("table_function_tensor"));
        auto tensor = std::static_pointer_cast< FieldConstant<3, FieldValue<3>::TensorFixed> >(tensor_base);
    	arma::mat::fixed<3,3> result;

    	result = tensor->time_dependent_value(0.2);
        EXPECT_DOUBLE_EQ( 1.4, result[0]);
        EXPECT_DOUBLE_EQ( 0.0, result[1]);
        EXPECT_DOUBLE_EQ( 1.2, result[2]);
        EXPECT_DOUBLE_EQ( 3.0, result[3]);
        EXPECT_DOUBLE_EQ( 2.8, result[4]);
        EXPECT_DOUBLE_EQ( 6.2, result[5]);
        EXPECT_DOUBLE_EQ( 4.4, result[6]);
        EXPECT_DOUBLE_EQ( 4.2, result[7]);
        EXPECT_DOUBLE_EQ( 5.6, result[8]);

    	result = tensor->time_dependent_value(1.5);
        EXPECT_DOUBLE_EQ( 4.0, result[0]);
        EXPECT_DOUBLE_EQ( 1.0, result[1]);
        EXPECT_DOUBLE_EQ( 3.5, result[2]);
        EXPECT_DOUBLE_EQ( 3.0, result[3]);
        EXPECT_DOUBLE_EQ( 2.5, result[4]);
        EXPECT_DOUBLE_EQ( 6.5, result[5]);
        EXPECT_DOUBLE_EQ( 5.0, result[6]);
        EXPECT_DOUBLE_EQ( 6.0, result[7]);
        EXPECT_DOUBLE_EQ( 3.5, result[8]);

    }

}


