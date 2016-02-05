/*
 * multi_field_test.cpp
 *
 *  Created on: Dec 11, 2012
 *      Author: jb
 */


#include <flow_gtest.hh>

#include <fields/multi_field.hh>
#include <input/type_base.hh>
#include <input/reader_to_storage.hh>

#include <input/type_output.hh>

#include <iostream>
using namespace std;

FLOW123D_FORCE_LINK_IN_PARENT(field_constant)

string field_constant_input = R"JSON(
{
    common={ 
        TYPE="FieldConstant", 
        value=[1, 2, 3]
    },
    transposed=[ 
        { TYPE="FieldConstant", value=1},
        { TYPE="FieldConstant", value=2},
        { TYPE="FieldConstant", value=3}
    ]
}
)JSON";

TEST(TransposeTo, field_constant) {
	MultiField<3, FieldValue<3>::Scalar> empty_mf;
	Input::Type::Record in_rec = Input::Type::Record("MultiFieldTest","")
	    .declare_key("common", empty_mf.get_multifield_input_type(), Input::Type::Default::obligatory(),"" )
	    .declare_key("transposed", empty_mf.get_multifield_input_type(), Input::Type::Default::obligatory(),"" )
	    .close();

	Input::ReaderToStorage json_reader(field_constant_input, in_rec, Input::FileFormat::format_JSON);
	Input::Record input = json_reader.get_root_interface<Input::Record>();

	Input::Array common = input.val<Input::Array>("common");
	Input::Array transposed = input.val<Input::Array>("transposed");

	EXPECT_EQ(common.size(), transposed.size());

	auto it_c = common.begin<Input::AbstractRecord>();
	for (auto it_t = transposed.begin<Input::AbstractRecord>(); it_t != transposed.end(); ++it_t, ++it_c) {
		Input::Record rec_t = (*it_t);
		Input::Record rec_c = (*it_c);
		EXPECT_DOUBLE_EQ( rec_t.val<double>("value"), rec_c.val<double>("value") );
	}
}
