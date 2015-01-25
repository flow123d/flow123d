/*
 * multi_field_test.cpp
 *
 *  Created on: Dec 11, 2012
 *      Author: jb
 */


#include <flow_gtest.hh>

#include <fields/multi_field.hh>
#include <input/type_base.hh>
#include <input/json_to_storage.hh>

#include <input/type_output.hh>

#include <iostream>
using namespace std;

string field_constant_input = R"JSON(
{
    common={ 
        TYPE="FieldConstant", 
        component_names=[ "field_1", "field_2", "field_3" ],
        common={ TYPE="FieldConstant", value=["1", "2", "3"]}
    },
    transposed={ 
        TYPE="FieldConstant", 
        component_names=[ "field_1", "field_2", "field_3" ],
        components=[ { TYPE="FieldConstant", value="1"}, { TYPE="FieldConstant", value="2"}, { TYPE="FieldConstant", value="3"} ]
    }
}
)JSON";

TEST(TransposeTo, field_constant) {
	MultiField<3, FieldValue<3>::Scalar> empty_mf;
	Input::Type::Record in_rec("MultiFieldTest","");
	in_rec.declare_key("common", empty_mf.get_input_type(), Input::Type::Default::obligatory(),"" );
	in_rec.declare_key("transposed", empty_mf.get_input_type(), Input::Type::Default::obligatory(),"" );
	in_rec.finish();

	std::cout << Input::Type::OutputText(&in_rec) << std::endl;
	std::cout << "-------------------------------------------------------------" << std::endl;
	std::cout << Input::Type::OutputJSONTemplate(&in_rec) << std::endl;

	//Input::JSONToStorage json_reader(field_constant_input, in_rec);
}
