/*
 * field_time_function_test.cpp
 *
 *  Created on: Dec 18, 2012
 *      Author: jb
 */



#define FEAL_OVERRIDE_ASSERTS

#include <flow_gtest.hh>
#include <memory>


#include "fields/field_constant.hh"
#include "fields/field_time_function.hh"
#include "fields/table_function.hh"
#include "input/input_type.hh"
#include "input/accessors.hh"
#include "input/reader_to_storage.hh"



FLOW123D_FORCE_LINK_IN_PARENT(field_time_function)

typedef FieldAlgorithmBase<3, FieldValue<3>::Tensor > TensorField;
typedef FieldAlgorithmBase<3, FieldValue<3>::Vector > VectorField;
typedef FieldAlgorithmBase<3, FieldValue<0>::Scalar > ScalarField;


string field_time_function_input = R"YAML(
table_function_scalar: !FieldTimeFunction
  time_function:
    - - 0.0
      - 0.5
    - - 1.0
      - 1.0
    - - 2.0
      - 3.0
table_function_vector: !FieldTimeFunction
  time_function:
    - - 0.0
      - [0.5, 1.5, 2.0]
    - - 1.0
      - [1.0, 1.5, 1.0]
    - - 2.0
      - [3.0, 1.5, 5.0]
table_function_tensor: !FieldTimeFunction
  time_function:
    - - 0.0
      - [ [1,3,4], [0,3,4], [1,6,6] ]
    - - 1.0
      - [ [3,3,6], [0,2,5], [2,7,4] ]
    - - 2.0
      - [ [5,3,4], [2,3,7], [5,6,3] ]
)YAML";


FieldAlgoBaseInitData get_field_init_data(std::string field_name, double max_limit = std::numeric_limits<double>::max() ) {
	return FieldAlgoBaseInitData(field_name, 3, UnitSI::dimensionless(), std::make_pair(0, max_limit),
			(FieldFlag::declare_input & FieldFlag::equation_input & FieldFlag::allow_output) );
}


TEST(FieldTableFunction, table_function) {

    // setup FilePath directories
    FilePath::set_io_dirs(".",UNIT_TESTS_SRC_DIR,"",".");

    Input::Type::Record rec_type = Input::Type::Record("FieldTableFunctionTest","")
        .declare_key("table_function_scalar", ScalarField::get_input_type_instance(), Input::Type::Default::obligatory(),"" )
        .declare_key("table_function_vector", VectorField::get_input_type_instance(), Input::Type::Default::obligatory(),"" )
        .declare_key("table_function_tensor", TensorField::get_input_type_instance(), Input::Type::Default::obligatory(),"" )
        .close();

    // read input string
    Input::ReaderToStorage reader( field_time_function_input, rec_type, Input::FileFormat::format_YAML );
    Input::Record in_rec=reader.get_root_interface<Input::Record>();

    Space<3>::Point point;
    point(0)=1.0; point(1)=2.0; point(2)=3.0;
    ElementAccessor<3> elm;
    auto conversion_ptr = std::make_shared<TimeUnitConversion>();

    {
        std::shared_ptr< ScalarField > scalar_base =
        		ScalarField::function_factory(in_rec.val<Input::AbstractRecord>("table_function_scalar"), get_field_init_data("table_function_scalar"));
        auto scalar = std::static_pointer_cast< FieldTimeFunction<3, FieldValue<0>::Scalar> >(scalar_base);

        scalar->set_time(0.0);
        EXPECT_DOUBLE_EQ( scalar->value(point, elm), 0.5 );
        scalar->set_time(0.2);
    	EXPECT_DOUBLE_EQ( scalar->value(point, elm), 0.6 );
        scalar->set_time(0.4);
    	EXPECT_DOUBLE_EQ( scalar->value(point, elm), 0.7 );
        scalar->set_time(0.6);
    	EXPECT_DOUBLE_EQ( scalar->value(point, elm), 0.8 );
        scalar->set_time(0.8);
    	EXPECT_DOUBLE_EQ( scalar->value(point, elm), 0.9 );
        scalar->set_time(1.0);
    	EXPECT_DOUBLE_EQ( scalar->value(point, elm), 1.0 );
        scalar->set_time(1.5);
    	EXPECT_DOUBLE_EQ( scalar->value(point, elm), 2.0 );
        scalar->set_time(2.0);
    	EXPECT_DOUBLE_EQ( scalar->value(point, elm), 3.0 );
        scalar->set_time(3.0);
    	EXPECT_DOUBLE_EQ( scalar->value(point, elm), 3.0 );
    }

    {
        std::shared_ptr< VectorField > vector_base =
        		VectorField::function_factory(in_rec.val<Input::AbstractRecord>("table_function_vector"), get_field_init_data("table_function_vector", 2.5));
        auto vector = std::static_pointer_cast< FieldTimeFunction<3, FieldValue<3>::Vector> >(vector_base);
    	arma::vec result;

    	vector->set_time(0.2);
    	result = vector->value(point, elm);
        EXPECT_DOUBLE_EQ( 0.6, result[0]);
        EXPECT_DOUBLE_EQ( 1.5, result[1]);
        EXPECT_DOUBLE_EQ( 1.8, result[2]);

    	vector->set_time(0.6);
    	result = vector->value(point, elm);
        EXPECT_DOUBLE_EQ( 0.8, result[0]);
        EXPECT_DOUBLE_EQ( 1.5, result[1]);
        EXPECT_DOUBLE_EQ( 1.4, result[2]);

    	vector->set_time(1.0);
    	result = vector->value(point, elm);
        EXPECT_DOUBLE_EQ( 1.0, result[0]);
        EXPECT_DOUBLE_EQ( 1.5, result[1]);
        EXPECT_DOUBLE_EQ( 1.0, result[2]);

    	vector->set_time(1.5);
    	result = vector->value(point, elm);
        EXPECT_DOUBLE_EQ( 2.0, result[0]);
        EXPECT_DOUBLE_EQ( 1.5, result[1]);
        EXPECT_DOUBLE_EQ( 3.0, result[2]);

    }

    {
        std::shared_ptr< TensorField > tensor_base =
        		TensorField::function_factory(in_rec.val<Input::AbstractRecord>("table_function_tensor"), get_field_init_data("table_function_tensor"));
        auto tensor = std::static_pointer_cast< FieldTimeFunction<3, FieldValue<3>::Tensor> >(tensor_base);
    	arma::mat::fixed<3,3> result;

    	tensor->set_time(0.2);
    	result = tensor->value(point, elm);
        EXPECT_DOUBLE_EQ( 1.4, result[0]);
        EXPECT_DOUBLE_EQ( 0.0, result[1]);
        EXPECT_DOUBLE_EQ( 1.2, result[2]);
        EXPECT_DOUBLE_EQ( 3.0, result[3]);
        EXPECT_DOUBLE_EQ( 2.8, result[4]);
        EXPECT_DOUBLE_EQ( 6.2, result[5]);
        EXPECT_DOUBLE_EQ( 4.4, result[6]);
        EXPECT_DOUBLE_EQ( 4.2, result[7]);
        EXPECT_DOUBLE_EQ( 5.6, result[8]);

    	tensor->set_time(1.5);
    	result = tensor->value(point, elm);
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


