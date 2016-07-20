/*
 * input_interface_test.cpp
 *
 *  Created on: Apr 4, 2012
 *      Author: jb
 *
 *
 */


#define FEAL_OVERRIDE_ASSERTS

#include <flow_gtest.hh>

#include "input/accessors.hh"
#include "input/input_type.hh"
#include "input/reader_to_storage.hh"

const string json_input = R"JSON(
{
	int_val = 1,
	big_val = 256,
	str_val = "some_string",
	color = "blue",
	time = [0.0, 10.0]
}
)JSON";


enum Colors {
    red = 0,
    green = 1,
    blue = 2
};


class InputTypeDispatchTest : public testing::Test, public Input::ReaderToStorage {
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

    	tpl = new Tuple("TimeTuple", "Definition of time.");
    	{
    		tpl->declare_key("start_time", Double(0.0), Default::obligatory(),
    				"Start of simulation");
    		tpl->declare_key("end_time", Double(1.0), Default::obligatory(),
    				"End of simulation");
    		tpl->close();
    	}

        root_record = new Record("RootRecord", "Root record of JSON input");
        {
        	root_record->declare_key("int_val", Integer(0), Default::obligatory(),
                 "Some integral value.");
        	root_record->declare_key("big_val", Integer(0), Default::obligatory(),
                 "Other integral value.");
        	root_record->declare_key("str_val", String(), Default::obligatory(),
                 "Some string value.");
        	root_record->declare_key("color", *sel, Default::obligatory(),
                 "Some enum value.");
        	root_record->declare_key("time", *tpl, Default::obligatory(),
                 "Some tuple value.");
        	root_record->close();
        }
    }
    virtual void TearDown() {
    };

    // overload parent class method in order to reset pointers
    void read_stream(istream &in, const Input::Type::TypeBase &root_type) {
    	this->storage_ = nullptr;
    	this->root_type_ = nullptr;
    	ReaderToStorage::read_stream(in, root_type, Input::FileFormat::format_JSON);
    }

    Input::Type::Record * root_record;
    Input::Type::Selection * sel;
    Input::Type::Tuple * tpl;
};


TEST_F(InputTypeDispatchTest, all) {

    Input::ReaderToStorage json_reader(json_input, *root_record, Input::FileFormat::format_JSON);
    Input::Record rec=json_reader.get_root_interface<Input::Record>();
    EXPECT_EQ(1, *(rec.find<int>("int_val")) );
    EXPECT_EQ("some_string", *(rec.find<std::string>("str_val")) );
    EXPECT_EQ(blue, *(rec.find<Colors>("color")) );
    Input::Iterator<Input::Tuple> it = rec.find<Input::Tuple>("time");
    EXPECT_FLOAT_EQ(0.0, *(it->find<double>("start_time")) );
    EXPECT_FLOAT_EQ(10.0, *(it->find<double>("end_time")) );

    EXPECT_THROW_WHAT({ rec.val<char>("big_val"); }, Input::ExcInputMessage, "Value out of bounds.");
}
