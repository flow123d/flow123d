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
#include <input/type_abstract.hh>
#include <input/type_selection.hh>
#include <input/reader_to_storage.hh>


using namespace Input;

class ReaderIncludeTest : public testing::Test, public Input::ReaderToStorage {
protected:

    virtual void SetUp() {
    }
    virtual void TearDown() {
    };

    static Type::Record & get_input_include_record() {
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

    static Type::Record & get_input_include_csv_record() {
    	static Type::Abstract abstr = Type::Abstract("Abstract", "Abstract type.")
				.close();

    	static Type::Record derived_rec = Type::Record("DerivedRec", "Some derived record.")
    			.derive_from( abstr )
    			.declare_key("some_double", Type::Double(), Type::Default::obligatory(), "desc.")
    			.declare_key("some_int", Type::Integer(), Type::Default::obligatory(), "desc.")
    			.close();

    	static Type::Record coords_rec = Type::Record("Coords", "3D coords.")
    			.declare_key("x", Type::Double(), Type::Default::obligatory(), "desc.")
    			.declare_key("y", Type::Double(), Type::Default::obligatory(), "desc.")
    			.declare_key("z", Type::Double(), Type::Default::obligatory(), "desc.")
    			.close();

    	static Type::Selection int_sel = Type::Selection("IntSelection")
    			.add_value(1,"one","")
				.add_value(2,"two","")
				.add_value(10,"ten","")
				.close();

    	static Type::Record sub_rec = Type::Record("SubRecord", "Sub record.")
    			.declare_key("coords", coords_rec, Type::Default::obligatory(), "desc.")
				.declare_key("abstract_key", abstr, Type::Default::obligatory(), "desc.")
				.declare_key("select", int_sel, Type::Default::obligatory(), "desc.")
    			.declare_key("a_key", Type::Integer(), Type::Default::obligatory(), "desc.")
    			.declare_key("b_key", Type::Double(), Type::Default::obligatory(), "desc.")
    			.declare_key("c_key", Type::Bool(), Type::Default::obligatory(), "desc.")
    			.declare_key("d_key", Type::String(), Type::Default::obligatory(), "desc.")
    			.declare_key("e_key", Type::Double(), Type::Default::obligatory(), "desc.")
    			.close();

    	return Type::Record("RootCsv", "Root record with CSV include.")
    			.declare_key("sub_rec", Type::Array(sub_rec), Type::Default::obligatory(), "desc.")
    			.declare_key("description", Type::String(), Type::Default::obligatory(), "desc.")
    			.close();
    };

    // overload parent class method in order to reset pointers
    void read_stream(istream &in, Type::TypeBase &root_type, FileFormat format = FileFormat::format_JSON) {
    	this->storage_ = nullptr;
    	this->root_type_ = nullptr;
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

    // check storage in test with including CSV input file
    void check_csv_storage() {
    	// reference values of integer, selection and string columns
    	static int ref_ints[] = {2,3,0,1,4,5};
    	static int ref_sels[] = {1,2,10,1,2,10};
    	static std::string ref_strs[] = {"ABC","DEF","GHI","JKL","MNO","PQR"};

    	EXPECT_NE((void *)NULL, storage_);
    	EXPECT_EQ(3, storage_->get_array_size());
    	EXPECT_EQ(6, storage_->get_item(1)->get_array_size());
    	for (unsigned int i=0; i<storage_->get_item(1)->get_array_size(); ++i) {
    		StorageBase *i_storage = storage_->get_item(1)->get_item(i);
    		EXPECT_EQ(9,                 i_storage->get_array_size());
    		EXPECT_EQ("SubRecord",       i_storage->get_item(0)->get_string());
    		EXPECT_EQ(4,                 i_storage->get_item(1)->get_array_size());
    		EXPECT_EQ("Coords",          i_storage->get_item(1)->get_item(0)->get_string());
    		EXPECT_DOUBLE_EQ(0.0,        i_storage->get_item(1)->get_item(1)->get_double());
    		EXPECT_DOUBLE_EQ(0.1*i,      i_storage->get_item(1)->get_item(2)->get_double());
    		EXPECT_DOUBLE_EQ(0.2*(i+1),  i_storage->get_item(1)->get_item(3)->get_double());
    		EXPECT_EQ(3,                 i_storage->get_item(2)->get_array_size());
    		EXPECT_EQ("DerivedRec",      i_storage->get_item(2)->get_item(0)->get_string());
    		EXPECT_DOUBLE_EQ(0.1*(i+1),  i_storage->get_item(2)->get_item(1)->get_double());
    		EXPECT_EQ(ref_ints[i]+3,     i_storage->get_item(2)->get_item(2)->get_int());
    		EXPECT_EQ(ref_sels[i],       i_storage->get_item(3)->get_int());
    		EXPECT_EQ(ref_ints[i],       i_storage->get_item(4)->get_int());
    		EXPECT_DOUBLE_EQ(0.5*i-0.5,  i_storage->get_item(5)->get_double());
    		EXPECT_FALSE(                i_storage->get_item(6)->get_bool());
    		EXPECT_EQ(ref_strs[i],       i_storage->get_item(7)->get_string());
    		EXPECT_DOUBLE_EQ(0.25*i-1.0, i_storage->get_item(8)->get_double());
    	}
    	EXPECT_EQ("RootCsv", storage_->get_item(0)->get_string() );
    	EXPECT_EQ("Example of CSV include", storage_->get_item(2)->get_string() );
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

const string import_yaml_auto_conversion = R"YAML(
sub_rec: !include
  include_simple.yaml
int_key: 10
)YAML";

const string import_json_auto_conversion = R"YAML(
sub_rec: !include
  include_simple.con
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
    read_stream(ss, ReaderIncludeTest::get_input_include_record(), FileFormat::format_YAML);
    check_simple_storage();
}

TEST_F(ReaderIncludeTest, include_json_to_yaml) {
    stringstream ss( import_json_to_yaml );
    read_stream(ss, ReaderIncludeTest::get_input_include_record(), FileFormat::format_YAML);
    check_simple_storage();
}

TEST_F(ReaderIncludeTest, include_yaml_auto_conversion) {
    stringstream ss( import_yaml_auto_conversion );
    read_stream(ss, ReaderIncludeTest::get_input_include_record(), FileFormat::format_YAML);
    check_simple_storage();
}

TEST_F(ReaderIncludeTest, include_json_auto_conversion) {
    stringstream ss( import_json_auto_conversion );
    read_stream(ss, ReaderIncludeTest::get_input_include_record(), FileFormat::format_YAML);
    check_simple_storage();
}

TEST_F(ReaderIncludeTest, include_yaml_to_json) {
    stringstream ss( import_yaml_to_json );
    read_stream(ss, ReaderIncludeTest::get_input_include_record(), FileFormat::format_JSON);
    check_simple_storage();
}

TEST_F(ReaderIncludeTest, include_json_to_json) {
    stringstream ss( import_json_to_json );
    read_stream(ss, ReaderIncludeTest::get_input_include_record(), FileFormat::format_JSON);
    check_simple_storage();
}


const string import_csv_to_yaml = R"YAML(
sub_rec: !include_csv
  file: include_simple.csv
  format:
    coords:
      x: 1
      y: 2
      z: 3
    abstract_key: !DerivedRec
      some_double: 9
      some_int: 10
    a_key: 0
    b_key: 4
    c_key: 6
    d_key: 7
    e_key: 5
    select: 11
description: Example of CSV include
)YAML";

const string import_csv_to_json = R"JSON(
{
  sub_rec = { 
    TYPE = "include_csv", 
    file = "include_simple.csv",
    format = { 
      coords={x=1, y=2, z=3}, 
      abstract_key={TYPE="DerivedRec", some_double=9, some_int=10}, 
      a_key=0, b_key=4, c_key=6, d_key=7, e_key=5, select=11
    }
  },
  description = "Example of CSV include"
}
)JSON";

TEST_F(ReaderIncludeTest, include_csv_to_yaml) {
    stringstream ss( import_csv_to_yaml );
    read_stream(ss, ReaderIncludeTest::get_input_include_csv_record(), FileFormat::format_YAML);
    check_csv_storage();
}

TEST_F(ReaderIncludeTest, include_csv_to_json) {
    stringstream ss( import_csv_to_json );
    read_stream(ss, ReaderIncludeTest::get_input_include_csv_record(), FileFormat::format_JSON);
    check_csv_storage();
}
