#include <iostream>
#include <fstream>
#include <cassert>
#include <cstdlib>
#include <cstdio>
#include <cfloat>


#include <gtest/gtest.h>

#include "json_spirit.h"
#include "input.hpp" //jediny potrebny header, vse ostatni se includuje vevnitr
#include "fjson_data.hpp"

using namespace std;
using namespace flow;


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

