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

TEST(flow_json_parser_stream, trivial_pure) {
    //read from stream, data from string
    const string data("{  \"flow\"  :  \"OK\"  }");

    istringstream in_s(data);

    ASSERT_TRUE( in_s.good() );

    Data_tree * tree = new Data_tree(in_s);
    ASSERT_FALSE( tree->err_status );

    stringstream ss1, ss2;
    ss1 << *tree;
    ss2 << tree->get_head();

    ASSERT_STREQ("{\"flow\":\"OK\"}", ss1.str().c_str());
    ASSERT_STREQ("{\"flow\":\"OK\"}", ss2.str().c_str());

    delete tree;
}

TEST(flow_json_parser_stream, trivial_wrongdata) {
    //read from stream, data from string
    const string data("blahblah");

    istringstream in_s(data);

    ASSERT_TRUE( in_s.good() );

    Data_tree * tree = new Data_tree(in_s);
    ASSERT_TRUE( tree->err_status );

    delete tree;
}

TEST(flow_json_parser_stream, json_pure) {
    //read from stream, data from multiline string
    istringstream in_s(flow_mini_json);
    ASSERT_TRUE( in_s.good() );
    Data_tree * tree = new Data_tree(in_s);
    ASSERT_FALSE( tree->err_status );

    delete tree;
}

TEST(flow_json_parser_stream, comments) {
    //read from stream, data from multiline string
    istringstream in_s(flow_json_comment_parser);
    ASSERT_TRUE( in_s.good() );
    Data_tree * tree = new Data_tree(in_s);
    ASSERT_FALSE( tree->err_status );

    delete tree;
}

TEST(flow_json_parser_stream, colon_equal) {
    //read from stream, data from multiline string
    istringstream in_s(flow_json_colon_eq);
    ASSERT_TRUE( in_s.good() );
    Data_tree * tree = new Data_tree(in_s);
    ASSERT_FALSE( tree->err_status );

    delete tree;
}

TEST(flow_json_parser_stream, quotes) {
    //read from stream, data from multiline string
    istringstream in_s(flow_json_quotes);
    ASSERT_TRUE( in_s.good() );
    Data_tree * tree = new Data_tree(in_s);
    ASSERT_FALSE( tree->err_status );

    delete tree;
}

TEST(flow_json_parser_stream, whitespace_separator) {
    //read from stream, data from multiline string
    istringstream in_s(flow_json_whitespace_separator);
    ASSERT_TRUE( in_s.good() );
    Data_tree * tree = new Data_tree(in_s);
    ASSERT_FALSE( tree->err_status );

    delete tree;
}

/*
 * **********************************************************************************************************************
 */
TEST(flow_json_parser_string, trivial_pure) {
    //read from string
    const string data("{  \"flow\"  :  \"OK\"  }");

    Data_tree * tree = new Data_tree(data);
    ASSERT_FALSE( tree->err_status );

    stringstream ss1, ss2;
    ss1 << *tree;
    ss2 << tree->get_head();

    ASSERT_STREQ("{\"flow\":\"OK\"}", ss1.str().c_str());
    ASSERT_STREQ("{\"flow\":\"OK\"}", ss2.str().c_str());

    delete tree;
}

TEST(flow_json_parser_string, trivial_wrongdata) {
    //read from string, data from string
    const string data("blahblah");

    Data_tree * tree = new Data_tree(data);
    ASSERT_TRUE( tree->err_status );

    delete tree;
}

TEST(flow_json_parser_string, json_pure) {
    //read from string, data from multiline string
    Data_tree * tree = new Data_tree(flow_mini_json);
    ASSERT_FALSE( tree->err_status );

    delete tree;
}

TEST(flow_json_parser_string, comments) {
    //read from string, data from multiline string
    Data_tree * tree = new Data_tree(flow_json_comment_parser);
    ASSERT_FALSE( tree->err_status );

    delete tree;
}

TEST(flow_json_parser_string, colon_equal) {
    //read from string, data from multiline string
    Data_tree * tree = new Data_tree(flow_json_colon_eq);
    ASSERT_FALSE( tree->err_status );

    delete tree;
}

TEST(flow_json_parser_string, quotes) {
    //read from string, data from multiline string
    Data_tree * tree = new Data_tree(flow_json_quotes);
    ASSERT_FALSE( tree->err_status );

    delete tree;
}

TEST(flow_json_parser_string, whitespace_separator) {
    //read from string, data from multiline string
    Data_tree * tree = new Data_tree(flow_json_whitespace_separator);
    ASSERT_FALSE( tree->err_status );

    delete tree;
}

/*
 * **********************************************************************************************************************
 */

TEST(Node_lib, access_to_non_existent) {

    //TODO
}

TEST(Node_lib, access_as_ancestor) {

    //TODO
}

/*
 * **********************************************************************************************************************
 */

int main(int argc, char **argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
