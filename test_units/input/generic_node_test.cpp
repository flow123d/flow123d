#include <iostream>
#include <fstream>
#include <cassert>
#include <cstdlib>
#include <cstdio>
#include <cfloat>


#include <gtest/gtest.h>

#include "json_spirit.h"
#include "input.hpp" //jediny potrebny header, vse ostatni se includuje vevnitr

TEST(json_spirit, all) {

}



using namespace std;
using namespace flow;

// Test data strings
string flow_mini_json = R"JSON(
{
    "flow_ini_version"  : "1.0",
    "comment"           : "verified by http://json.parser.online.fr/ to be VALID JSON",

    "global" : {
        "problem_type"  : 1,
        "description"   : "test1",
        "save_step"     : 0.1,
        "density_on"    : false,
        "nothing"       : null
    },

    "input" : {
        "file_type"         : 1,
        "mesh"              : "./input/test1.msh"
    },

    "constants" : {
        "g"       :9.81,
        "rho"     :1000
    },

    "sp" : {
        "drfl"              :1e-009,
        "l_size"            :80
    },

    "z_test_weird_array" : [ [0], { "a" : 1 }, 2, {}, [] ]
}
)JSON";

string flow_json_comma_spearator = R"JSON(
{
    "key1" = "value1",
    "key2" : "value2",
     key3  = "value3",
     key4  : "value4"
    "key5" = "value5"
    "key6" : "value6"
     key7  = "value7"
     key8  : "value8"
}
)JSON";



TEST(json_spirit_test, read) {

    json_spirit::mValue tree_root;

    //read from stream
    ifstream in_s(flow_mini_json);
    json_spirit::read(in_s, tree_root);
    in_s.close();

    // test root node type
    EXPECT_EQ( json_spirit::obj_type, tree_root.type() );

    //get root record - flow ini has as root Object, no need to test
    json_spirit::mObject &root = tree_root.get_obj();

    //with iterator helper
    json_spirit::mObject::iterator i = root.find("flow_ini_version");
    EXPECT_EQ( string("1.0"), i->second.get_str() );

    //template request
    EXPECT_EQ( string("1.0"), i->second.get_value<string>() );

    //with built-in
    EXPECT_EQ( "verified by http://json.parser.online.fr/ to be VALID JSON",
            root.find("comment")->second.get_str() );

    //deeper in hierarchy, other type
    EXPECT_EQ( 1, root.find("global")->second.get_obj().find("problem_type")->second.get_int() );
}

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

/*
bool minitest_stream( const string & fname );
bool minitest_string( const string & fname );

void minitest( const string & fname )
{
    cout << "string test: ";
    if ( minitest_string( fname ) )
        cout << "FAILED" << endl;
    else
        cout << "OK" << endl;


    cout << "stream test: ";
    if ( minitest_stream( fname ) )
        cout << "FAILED" << endl;
    else
        cout << "OK" << endl;
}

int main()
{




    //test JSON FLOW extended format
    cout << endl << "===== JSON FLOW extended format =====" << endl;

    cout << "Test filtru komentaru:" << endl;
    minitest("src/test/comments.fjson");

    cout << "Test carek:" << endl;
    minitest("src/test/carky.fjson");

    cout << "Test rovnitek:" << endl;
    minitest("src/test/dvojtecka-rovnitko.fjson");

    cout << "Test uvozovek:" << endl;
    minitest("src/test/uvozovky.fjson");

    cout << "END." << endl << flush;
    return 0;
}

bool minitest_stream( const string & fname )
{
    Data_tree * tree;
    string in_str;
    bool retval;

    ifstream in_s( fname.c_str() );
    if (in_s.is_open())
    {
        tree = new Data_tree(in_s);
        in_s.close();
    } else {
        cout << "Unable to open file!" << endl;
        return false;
    }

    retval = tree->err_status;

    delete tree;
    return retval;
}

bool minitest_string( const string & fname )
{
    Data_tree * tree;
    string in_str;
    bool retval;

    ifstream in_s( fname.c_str() );
    if (in_s.is_open())
    {
        while (in_s.good())
        {
            in_str.push_back(in_s.get());
        }
        in_s.close();
        // The complete file content is in memory now...
    } else {
        cout << "Unable to open file!" << endl;
        return false;
    }

    tree = new Data_tree(in_str);

    retval = tree->err_status;

    delete tree;
    return retval;
}
*/
