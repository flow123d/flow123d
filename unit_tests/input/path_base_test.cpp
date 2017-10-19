/*
 * path_base_test.cpp
 *
 *  Created on: May 7, 2012
 *      Author: jb
 */

/**
 * TODO: test catching of errors in JSON file format.
 */

#define FEAL_OVERRIDE_ASSERTS

#include <flow_gtest.hh>
#include <fstream>


#include "input/reader_to_storage.hh"
#include "input/path_base.hh"
#include "input/path_json.hh"
#include "input/path_yaml.hh"

using namespace std;

using namespace Input;


TEST(PathJSON, all) {
::testing::FLAGS_gtest_death_test_style = "threadsafe";

    ifstream in_str((string(UNIT_TESTS_SRC_DIR) + "/input/reader_to_storage_test.con").c_str());
    PathJSON path(in_str);

    { ostringstream os;
    os << path;
    EXPECT_EQ("/",os.str());
    }

    path.down(6);
    { ostringstream os;
    os << path;
    EXPECT_EQ("/6",os.str());
    }

    path.down("a");
    { ostringstream os;
    os << path;
    EXPECT_EQ("/6/a",os.str());
    }
    EXPECT_EQ(1,path.find_ref_node()->get_int_value() );

    path.go_to_root();
    path.down(6);
    path.down("b");
    EXPECT_EQ("ctyri",path.find_ref_node()->get_string_value() );
}


TEST(PathJSON, errors) {
::testing::FLAGS_gtest_death_test_style = "threadsafe";
	ifstream in_str((string(UNIT_TESTS_SRC_DIR) + "/input/reader_to_storage_test.con").c_str());

    PathJSON path(in_str);
    path.down(8);
    EXPECT_THROW_WHAT( { path.find_ref_node();}, PathBase::ExcRefOfWrongType,"has wrong type, should by string." );

    path.go_to_root();
    path.down(9); // "REF":"/5/10"
    EXPECT_THROW_WHAT( { path.find_ref_node();}, PathBase::ExcReferenceNotFound, "index out of size of Array" );

    path.go_to_root();
    path.down(10); // "REF":"/6/../.."
    EXPECT_THROW_WHAT( { path.find_ref_node();}, PathBase::ExcReferenceNotFound, "can not go up from root" );

    path.go_to_root();
    path.down(11); // "REF":"/key"
    EXPECT_THROW_WHAT( { path.find_ref_node();}, PathBase::ExcReferenceNotFound, "there should be Record" );

    path.go_to_root();
    path.down(12); // "REF":"/6/key"
    EXPECT_THROW_WHAT( { path.find_ref_node();}, PathBase::ExcReferenceNotFound, "key 'key' not found" );
}



TEST(PathYAML, all) {
::testing::FLAGS_gtest_death_test_style = "threadsafe";

	ifstream in_str((string(UNIT_TESTS_SRC_DIR) + "/input/reader_to_storage_test.yaml").c_str());
	PathYAML path(in_str);

    { ostringstream os;
    os << path;
    EXPECT_EQ("/",os.str());
    }

    path.down(6);
    { ostringstream os;
    os << path;
    EXPECT_EQ("/6",os.str());
    }

    path.down("a");
    { ostringstream os;
    os << path;
    EXPECT_EQ("/6/a",os.str());
    }

    path.up();
    { ostringstream os;
    os << path;
    EXPECT_EQ("/6",os.str());
    }

    path.down("b");
    { ostringstream os;
    os << path;
    EXPECT_EQ("/6/b",os.str());
    }

}


TEST(PathYAML, values) {
::testing::FLAGS_gtest_death_test_style = "threadsafe";

	ifstream in_str((string(UNIT_TESTS_SRC_DIR) + "/input/reader_to_storage_test.yaml").c_str());
	PathYAML path(in_str);

	path.down(0); // bool value
	EXPECT_FALSE(path.get_bool_value());
	EXPECT_THROW( { path.get_int_value(); }, ReaderToStorage::ExcInputError );
	EXPECT_THROW( { path.get_double_value(); }, ReaderToStorage::ExcInputError );
	EXPECT_STREQ("false", path.get_string_value().c_str());
	path.up();

	path.down(1); // int value
	EXPECT_THROW( { path.get_bool_value(); }, ReaderToStorage::ExcInputError );
	EXPECT_EQ(1, path.get_int_value());
	EXPECT_FLOAT_EQ(1.0, path.get_double_value());
	EXPECT_STREQ("1", path.get_string_value().c_str());
	path.up();

	path.down(3); // double value
	EXPECT_THROW( { path.get_bool_value(); }, ReaderToStorage::ExcInputError );
	EXPECT_THROW( { path.get_int_value(); }, ReaderToStorage::ExcInputError );
	EXPECT_FLOAT_EQ(3.3, path.get_double_value());
	EXPECT_STREQ("3.3", path.get_string_value().c_str());
	path.up();

	path.down(4); // string value
	EXPECT_THROW( { path.get_bool_value(); }, ReaderToStorage::ExcInputError );
	EXPECT_THROW( { path.get_int_value(); }, ReaderToStorage::ExcInputError );
	EXPECT_THROW( { path.get_double_value(); }, ReaderToStorage::ExcInputError );
	EXPECT_STREQ("ctyri", path.get_string_value().c_str());
	path.up();

	path.down(9); // int64 value
	EXPECT_EQ(5000000000000, path.get_int_value());
	path.up();

	path.down(6); // record
	std::set<std::string> set;
	path.get_record_key_set(set);
	EXPECT_EQ(2, set.size());
	EXPECT_TRUE( set.find("a")!=set.end() );
	EXPECT_TRUE( set.find("b")!=set.end() );
	EXPECT_FALSE( set.find("c")!=set.end() );

	path.down("b"); // reference
	EXPECT_STREQ("ctyri", path.get_string_value().c_str());
}
