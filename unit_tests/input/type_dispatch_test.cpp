/*
 * input_interface_test.cpp
 *
 *  Created on: Apr 4, 2012
 *      Author: jb
 *
 *
 */


#include <flow_gtest.hh>

#include "input/accessors.hh"
#include "input/input_type.hh"
#include "input/json_to_storage.hh"

const string json_input = R"JSON(
{
	int_val = 1,
	str_val = "some_string",
	color = "blue"
}
)JSON";


enum Colors {
    red = 0,
    green = 1,
    blue = 2
};


class InputTypeDispatchTest : public testing::Test, public Input::JSONToStorage {
protected:

    virtual void SetUp() {
    	using namespace Input::Type;

    	sel = new Selection("Colors");
    	{
    	    sel->add_value(red, "red");
    	    sel->add_value(green, "green");
    	    sel->add_value(blue, "blue");
    	    sel->close();
    	}

        root_record = new Record("RootRecord", "Root record of JSON input");
        {
        	root_record->declare_key("int_val", Integer(0), Default::obligatory(),
                 "Some integral value.");
        	root_record->declare_key("str_val", String(), Default::obligatory(),
                 "Some string value.");
        	root_record->declare_key("color", *sel, Default::obligatory(),
                 "Some enum value.");
        	root_record->close();
        }
    }
    virtual void TearDown() {
    };

    // overload parent class method in order to reset pointers
    void read_stream(istream &in, const Input::Type::TypeBase &root_type) {
    	this->storage_ = nullptr;
    	this->root_type_ = nullptr;
    	JSONToStorage::read_stream(in,root_type);
    }

    Input::Type::Record * root_record;
    Input::Type::Selection * sel;
};


TEST_F(InputTypeDispatchTest, all) {

    Input::JSONToStorage json_reader(json_input, *root_record);
    Input::Record rec=json_reader.get_root_interface<Input::Record>();
    EXPECT_EQ(1, *(rec.find<int>("int_val")) );
    EXPECT_EQ("some_string", *(rec.find<std::string>("str_val")) );
    //EXPECT_EQ("blue", *(rec.find<Input::Enum>("color")) );
}
