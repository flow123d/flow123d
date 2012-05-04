/*
 * input_reader_test.cpp
 *
 *  Created on: 25.4.2012
 *      Author: jiri
 */
#include <iostream>
#include <sstream>
#include <fstream>
#include <string>

#include <gtest/gtest.h>

#include "json_spirit.h"
#include "input.hpp"
#include "fjson_data.hpp"

using namespace std;
using namespace flow;

/*

Fatal assertion                 Nonfatal assertion          Verifies
ASSERT_TRUE(condition);         EXPECT_TRUE(condition);     condition is true
ASSERT_FALSE(condition);        EXPECT_FALSE(condition);    condition is false

ASSERT_EQ(expected, actual);    EXPECT_EQ(expected, actual);    expected == actual
ASSERT_NE(val1, val2);          EXPECT_NE(val1, val2);          val1 != val2
ASSERT_LT(val1, val2);          EXPECT_LT(val1, val2);          val1 < val2
ASSERT_LE(val1, val2);          EXPECT_LE(val1, val2);          val1 <= val2
ASSERT_GT(val1, val2);          EXPECT_GT(val1, val2);          val1 > val2
ASSERT_GE(val1, val2);          EXPECT_GE(val1, val2);          val1 >= val2

ASSERT_STREQ(expected_str, actual_str);     EXPECT_STREQ(expected_str, actual_str);     the two C strings have the same content
ASSERT_STRNE(str1, str2);                   EXPECT_STRNE(str1, str2);                   the two C strings have different content
ASSERT_STRCASEEQ(expected_str, actual_str); EXPECT_STRCASEEQ(expected_str, actual_str); the two C strings have the same content, ignoring case
ASSERT_STRCASENE(str1, str2);               EXPECT_STRCASENE(str1, str2);               the two C strings have different content, ignoring case

EXPECT_DEATH()
ASSERT_DEATH()
 */

/*
 * **********************************************************************************************************************
 */

TEST(json_spirit_test, read) {
    json_spirit::mValue tree_root;

    //read from stream, data from multiline string
    {
        istringstream in_s(flow_mini_json);
        ASSERT_TRUE( in_s.good() );
        ASSERT_TRUE( json_spirit::read(in_s, tree_root) );
    }

    // test root node type
    ASSERT_EQ( json_spirit::obj_type, tree_root.type() );

    //get root record - flow ini has as root Object, no need to test
    json_spirit::mObject &root = tree_root.get_obj();

    //with iterator helper
    json_spirit::mObject::iterator i = root.find("flow_ini_version");
    ASSERT_STREQ( "1.0", i->second.get_str().c_str() );

    //template request
    ASSERT_STREQ( "1.0", i->second.get_value<string>().c_str() );

    //with built-in
    ASSERT_STREQ( "verified by http://json.parser.online.fr/ to be VALID JSON", root.find("comment")->second.get_str().c_str() );

    //with built-in, template
    ASSERT_STREQ( "verified by http://json.parser.online.fr/ to be VALID JSON", root.find("comment")->second.get_value<string>().c_str() );


    //deeper in hierarchy, other type
    ASSERT_EQ( 1, root.find("global")->second.get_obj().find("problem_type")->second.get_int() );
}

/*
 * **********************************************************************************************************************
 */

TEST(Node_lib, access_to_non_existent) {
    //non-existent node example
    const int default_int = 12321;

    Record_node rnode;
    Generic_node & gnode_r = rnode;

    //if does not exist, return dafault
    ASSERT_EQ( default_int, gnode_r.get_key("foo").get_item(10).as_value().get_int(default_int) );


    //if does not exist, return zero and errcode
    int errcode;
    ASSERT_EQ( 0, gnode_r.get_key("foo").get_item(10).as_value().get_int_check(errcode) );
    ASSERT_NE( 0, errcode );

    //if does not exist, die
    ASSERT_DEATH( gnode_r.get_key("foo").get_item(10).as_value().get_int(), "Internal Error" );
}

TEST(Node_lib, access_as_ancestor) {
    //access to instance as an ancestor
    const int value(123);

    Value_node vnode(value);
    Generic_node & gnode_r = vnode;
    Generic_node * gnode_p = &vnode;

    ASSERT_EQ( value, vnode.get_int() );
    ASSERT_EQ( value, gnode_r.as_value().get_int() );
    ASSERT_EQ( value, (*gnode_p).as_value().get_int() );
}

TEST(Node_lib, copy_node_value) {
    //try all to copy all node types, check independency
    //value
    //set first
    const int int_def(123);
    Value_node vnode(int_def);
    ASSERT_EQ( type_number, vnode.get_type() );
    ASSERT_EQ( int_def, vnode.get_int() );

    //copy second
    Value_node vnode2 = vnode;
    ASSERT_EQ( type_number, vnode2.get_type() );
    ASSERT_EQ( int_def, vnode2.get_int() );

    //change first
    const int int_def2(555);
    const bool bool_def(true);
    vnode.set_value(int_def2); //first change value
    vnode.set_value(bool_def); //and than both value and type
    ASSERT_EQ( type_bool, vnode.get_type() );
    ASSERT_EQ( bool_def, vnode.get_bool() );
    ASSERT_EQ( type_number, vnode2.get_type() );
    ASSERT_EQ( int_def, vnode2.get_int() );
}

TEST(Node_lib, copy_node_record) {
    //record - create record with one int
    //set first
    const int int_def(123);
    Value_node vnode(int_def);
    Record_node rnode;
    ASSERT_EQ( type_number, vnode.get_type() );
    ASSERT_EQ( int_def, vnode.get_int() );
    ASSERT_EQ( type_record, rnode.get_type() );
    ASSERT_EQ( 0, rnode.get_record_size() );

    const string key_name("int");
    rnode.insert_key( key_name , &vnode );
    ASSERT_EQ( 1, rnode.get_record_size() );
    ASSERT_EQ( type_number, rnode.get_key(key_name).get_type() );
    ASSERT_EQ( int_def, rnode.get_key(key_name).as_value().get_int() );

    //copy second
    Record_node rnode2 = rnode;
    ASSERT_EQ( 1, rnode2.get_record_size() );
    ASSERT_EQ( type_number, rnode2.get_key(key_name).get_type() );
    ASSERT_EQ( int_def, rnode2.get_key(key_name).as_value().get_int() );

    //change first
    const int int_def2(555);
    const bool bool_def(true);
    rnode.get_key(key_name).as_value().set_value(int_def2); //first change value
    rnode.get_key(key_name).as_value().set_value(bool_def); //and than both value and type
    ASSERT_EQ( type_bool, rnode.get_key(key_name).get_type() );
    ASSERT_EQ( bool_def, rnode.get_key(key_name).as_value().get_bool() );
    ASSERT_EQ( type_number, rnode2.get_key(key_name).get_type() );
    ASSERT_EQ( int_def, rnode2.get_key(key_name).as_value().get_int() );
}

TEST(Node_lib, copy_node_vector) {
    //vector - create vector with one int
    //set first
    const int int_def(123);
    Value_node vnode(int_def);
    Vector_node vecnode;
    ASSERT_EQ( type_number, vnode.get_type() );
    ASSERT_EQ( int_def, vnode.get_int() );
    ASSERT_EQ( type_vector, vecnode.get_type() );
    ASSERT_EQ( 0, vecnode.get_array_size() );

//    vecnode.insert_item( 0 , vnode );
//    ASSERT_EQ( 1, vecnode.get_array_size() );
//    ASSERT_EQ( type_number, vecnode.get_item(0).get_type() );
//    ASSERT_EQ( int_def, vecnode.get_item(0).as_value().get_int() );
//
//    //copy second
//    Vector_node vecnode2 = vecnode;
//    ASSERT_EQ( 1, vecnode2.get_array_size() );
//    ASSERT_EQ( type_number, vecnode2.get_item(0).get_type() );
//    cout << "1" << endl ;
//    ASSERT_EQ( int_def, vecnode2.get_item(0).as_value().get_int() );
//    cout << "2";
//
//    //change first
//    const int int_def2(555);
//    const bool bool_def(true);
//    vecnode.get_item(0).as_value().set_value(int_def2); //first change value
//    vecnode.get_item(0).as_value().set_value(bool_def); //and than both value and type
//    ASSERT_EQ( type_bool, vecnode.get_item(0).get_type()  );
//    ASSERT_EQ( bool_def, vecnode.get_item(0).as_value().get_bool() );
//    ASSERT_EQ( type_number, vecnode2.get_item(0).get_type()  );
//    ASSERT_EQ( int_def, vecnode2.get_item(0).as_value().get_int() );
}

/*
 * **********************************************************************************************************************
 */

TEST(flow_json_parser_stream, trivial_pure) {
    //read from stream, data from string
    const string data("{  \"flow\"  :  \"OK\"  }");

    istringstream in_s(data);

    ASSERT_TRUE( in_s.good() );

    Data_tree tree(in_s);
    ASSERT_FALSE( tree.err_status );

    stringstream ss1, ss2;
    ss1 << tree;
    ss2 << tree.get_head();

    ASSERT_STREQ("{\"flow\":\"OK\"}", ss1.str().c_str());
    ASSERT_STREQ("{\"flow\":\"OK\"}", ss2.str().c_str());
}

TEST(flow_json_parser_stream, trivial_wrongdata) {
    //read from stream, data from string
    const string data("blahblah");

    istringstream in_s(data);
    ASSERT_TRUE( in_s.good() );

    Data_tree tree(in_s);
    ASSERT_TRUE( tree.err_status );
}

TEST(flow_json_parser_stream, json_pure) {
    //read from stream, data from multiline string
    istringstream in_s(flow_mini_json);
    ASSERT_TRUE( in_s.good() );
    Data_tree tree(in_s);
    ASSERT_FALSE( tree.err_status );
}

TEST(flow_json_parser_stream, comments) {
    //read from stream, data from multiline string
    istringstream in_s(flow_json_comment_parser);
    ASSERT_TRUE( in_s.good() );
    Data_tree tree(in_s);
    ASSERT_FALSE( tree.err_status );
}

TEST(flow_json_parser_stream, colon_equal) {
    //read from stream, data from multiline string
    istringstream in_s(flow_json_colon_eq);
    ASSERT_TRUE( in_s.good() );
    Data_tree tree(in_s);
    ASSERT_FALSE( tree.err_status );
}

TEST(flow_json_parser_stream, quotes) {
    //read from stream, data from multiline string
    istringstream in_s(flow_json_quotes);
    ASSERT_TRUE( in_s.good() );
    Data_tree tree(in_s);
    ASSERT_FALSE( tree.err_status );
}

TEST(flow_json_parser_stream, whitespace_separator) {
    //read from stream, data from multiline string
    istringstream in_s(flow_json_whitespace_separator);
    ASSERT_TRUE( in_s.good() );
    Data_tree tree(in_s);
    ASSERT_FALSE( tree.err_status );
}

/*
 * **********************************************************************************************************************
 */

TEST(flow_json_parser_string, trivial_pure) {
    //read from string
    const string data("{  \"flow\"  :  \"OK\"  }");

    Data_tree tree(data);
    ASSERT_FALSE( tree.err_status );

    stringstream ss1, ss2;
    ss1 << tree;
    ss2 << tree.get_head();

    ASSERT_STREQ("{\"flow\":\"OK\"}", ss1.str().c_str());
    ASSERT_STREQ("{\"flow\":\"OK\"}", ss2.str().c_str());
}

TEST(flow_json_parser_string, trivial_wrongdata) {
    //read from string, data from string
    Data_tree tree("blahblah");
    ASSERT_TRUE( tree.err_status );
}

TEST(flow_json_parser_string, json_pure) {
    //read from string, data from multiline string
    Data_tree tree(flow_mini_json);
    ASSERT_FALSE( tree.err_status );
}

TEST(flow_json_parser_string, comments) {
    //read from string, data from multiline string
    Data_tree tree(flow_json_comment_parser);
    ASSERT_FALSE( tree.err_status );
}

TEST(flow_json_parser_string, colon_equal) {
    //read from string, data from multiline string
    Data_tree tree(flow_json_colon_eq);
    ASSERT_FALSE( tree.err_status );
}

TEST(flow_json_parser_string, quotes) {
    //read from string, data from multiline string
    Data_tree tree(flow_json_quotes);
    ASSERT_FALSE( tree.err_status );
}

TEST(flow_json_parser_string, whitespace_separator) {
    //read from string, data from multiline string
    Data_tree tree(flow_json_whitespace_separator);
    ASSERT_FALSE( tree.err_status );
}

/*
 * **********************************************************************************************************************
 */
/*
TEST(data_tree_test, read) {
    ifstream in_s(flow_mini_json);
    Data_tree flow_tree(in_s);
    in_s.close();

    // check correct read
    EXPECT_EQ( 0, flow_tree.err_status);

    // check existing key
    Generic_node & root_node = flow_tree.get_head();
    EXPECT_EQ( "1.0", root_node.get_key("flow_ini_version").get_string());
    // non-existing key with default value
    EXPECT_EQ( "0.0", root_node.get_key("flow_ini_verssion").get_string("0.0"));

    // non-existent key with check
    int err;
    root_node.get_key("flow_ini_verssion").get_string_check(err);
    EXPECT_NE(0, err);

    int default_int = 12321;
    int returned_int;
    returned_int = root_node.get_key("foo").get_item(10).as_value().get_int(default_int);
    EXPECT_EQ(default_int, returned_int);

    // report error for nonexistent key
    EXPECT_DEATH(root_node.get_key("flow_ini_verssion").get_string(), "some error message\n");


}
*/
/*
 * **********************************************************************************************************************
 */

int main(int argc, char **argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
