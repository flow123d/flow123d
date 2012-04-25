/*
 * input_reader_test.cpp
 *
 *  Created on: 25.4.2012
 *      Author: jiri
 */
#include <iostream>
#include <fstream>

#include <gtest/gtest.h>

#include "input/input.hpp"

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


TEST(JSON_spirit,read_stream) {
    //using namespace std;
    json_spirit::mValue tree_root;
    ifstream in_s("flow_mini.json");
    ASSERT_FALSE( in_s.fail() );
    ASSERT_TRUE( json_spirit::read(in_s, tree_root) );
    in_s.close();
    ASSERT_FALSE( in_s.fail() );
}

int main(int argc, char **argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
