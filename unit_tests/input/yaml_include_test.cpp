/*
 * yaml_include_test.cpp
 *
 *  Created on: Apr 4, 2012
 *      Author: jb
 *

 *
 */


#define FEAL_OVERRIDE_ASSERTS

#include <flow_gtest.hh>
#include <input/type_record.hh>
#include <input/reader_to_storage.hh>


using namespace Input;

class ReaderIncludeTest : public testing::Test, public Input::ReaderToStorage {
protected:

    virtual void SetUp() {
    }
    virtual void TearDown() {
    };

    static Type::Record & get_input_record() {
    	static Type::Record sub_rec = Type::Record("SubRecord", "Sub record.")
    			.declare_key("a_key", Type::Array( Type::Double() ), Type::Default::obligatory(), "desc.")
    			.declare_key("b_key", Type::Array( Type::Integer() ), Type::Default::obligatory(), "desc.")
    			.declare_key("c_key", Type::Bool(), Type::Default::obligatory(), "desc.")
    			.close();

    	return Type::Record("Root", "Root record.")
    			.declare_key("sub_rec", sub_rec, Type::Default::obligatory(), "desc.")
    			.declare_key("int_key", Type::Integer(), Type::Default::obligatory(), "desc.")
    			.close();
    };

    // overload parent class method in order to reset pointers
    void read_stream(istream &in, Type::TypeBase &root_type, FileFormat format = FileFormat::format_JSON) {
    	this->storage_ = nullptr;
    	this->root_type_ = nullptr;
    	this->try_transpose_read_ = false;
    	root_type.finish();
    	ReaderToStorage::read_stream(in, root_type, format);
    };

    // check storage in test with including YAML or JSON input file
    void check_simple_storage() {
        EXPECT_NE((void *)NULL, storage_);
        EXPECT_EQ(3,   storage_->get_array_size());
        EXPECT_EQ(4,   storage_->get_item(1)->get_array_size() );
        EXPECT_EQ(3,   storage_->get_item(1)->get_item(1)->get_array_size() );
        EXPECT_EQ(1.0, storage_->get_item(1)->get_item(1)->get_item(0)->get_double() );
        EXPECT_EQ(3,   storage_->get_item(1)->get_item(2)->get_array_size() );
        EXPECT_EQ(5,   storage_->get_item(1)->get_item(2)->get_item(1)->get_int() );
        EXPECT_TRUE(   storage_->get_item(1)->get_item(3)->get_bool() );
        EXPECT_EQ(10,  storage_->get_item(2)->get_int() );
    };
};


const string import_yaml_to_yaml = R"YAML(
sub_rec: !include
  file: include_simple.yaml
int_key: 10
)YAML";

const string import_json_to_yaml = R"YAML(
sub_rec: !include
  file: include_simple.con
int_key: 10
)YAML";

const string import_yaml_to_json = R"JSON(
{
  sub_rec = { TYPE="include", file="include_simple.yaml" },
  int_key = 10
}
)JSON";

const string import_json_to_json = R"JSON(
{
  sub_rec = { TYPE="include", file="include_simple.con" },
  int_key = 10
}
)JSON";


TEST_F(ReaderIncludeTest, include_yaml_to_yaml) {
    stringstream ss( import_yaml_to_yaml );
    read_stream(ss, ReaderIncludeTest::get_input_record(), FileFormat::format_YAML);
    check_simple_storage();
}

TEST_F(ReaderIncludeTest, include_json_to_yaml) {
    stringstream ss( import_json_to_yaml );
    read_stream(ss, ReaderIncludeTest::get_input_record(), FileFormat::format_YAML);
    check_simple_storage();
}

TEST_F(ReaderIncludeTest, include_yaml_to_json) {
    stringstream ss( import_yaml_to_json );
    read_stream(ss, ReaderIncludeTest::get_input_record(), FileFormat::format_JSON);
    check_simple_storage();
}

TEST_F(ReaderIncludeTest, include_json_to_json) {
    stringstream ss( import_json_to_json );
    read_stream(ss, ReaderIncludeTest::get_input_record(), FileFormat::format_JSON);
    check_simple_storage();
}
