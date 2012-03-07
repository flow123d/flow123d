#include <iostream>
#include <fstream>
#include <cassert>
#include <cstdlib>
#include <cstdio>
#include <cfloat>

#include "../json_spirit/json_spirit.h"
#include "../input.hpp" //jediny potrebny header, vse ostatni se includuje vevnitr

using namespace std;
using namespace flow;


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
    cout << "===== JSON direct =====" << endl;

    //read from stream
    json_spirit::mValue tree_root;
    ifstream in_s("src/test/flow_mini.json");
    json_spirit::read(in_s, tree_root);
    in_s.close();

    //print all
    cout << json_spirit::write(tree_root) << endl << endl;

    //root record, or array??
    switch ( tree_root.type() ) {
    case json_spirit::obj_type:
        cout << "nalezen Record" << endl;
        break;
    case json_spirit::array_type:
        cout << "nalezeno Array" << endl;
        break;
    default:
        cout << "wtf?" << endl;
        break;
    }

    //print just some
    //get root record - flow ini has as root Object, no need to test
    //(the other root type can be array)
    json_spirit::mObject &root = tree_root.get_obj();

    //with iterator helper
    json_spirit::mObject::iterator i;
    i = root.find("flow_ini_version");
    cout << "verze: \"" << i->second.get_str() << "\"" << endl;

    //template request
    cout << "verze: \"" << i->second.get_value<string>() << "\"" << endl;

    //with built-in
    cout << "komentar: \"" << root.find("comment")->second.get_str() << "\"" << endl;
    cout << "komentar: \"" << root.find("comment")->second.get_value<string>() << "\"" << endl;

    //deeper in hierarchy
    cout << "ini-global-description: \"" << root.find("global")->second.get_obj().find("description")->second.get_str() << "\"" << endl;
    cout << "ini-global-description: \"" << root.find("global")->second.get_obj().find("description")->second.get_value<string>() << "\"" << endl;

    //test Node library
    cout << endl << "===== NODE LIBRARY =====" << endl;
    {
        Data_tree * tree;

        //read from file
        ifstream in_s("src/test/flow_mini.json");
        tree = new Data_tree(in_s);
        in_s.close();

        cout << "Tree error: " << tree->err_status << endl;

        cout << "JSON dump: ";
        tree->tree_dump_json();
        cout << endl;

        cout << "JSON tree << dump :" << (*tree) << endl;

        Generic_node & nodes = tree->get_head();

        cout << "Node tree << dump:" << endl << nodes << endl;

        //non-existent node example
/*
        Generic_node gnode;
        cout << endl << "Default demo:" << endl;
        int default_int = 12321;
        int returned_int;
        returned_int = gnode.get_key("foo").get_item(10).as_value().get_int(default_int);
        cout << "default_int=" << default_int << endl;
        cout << "returned_int=" << returned_int << endl;
*/
        //access to instance as an ancestor
/*
        Value_node vnode;
        Generic_node & gnode_r = vnode;
        Generic_node * gnode_p;
        gnode_p = new Value_node;
        cout << endl << "Access as ancestor demo:" << endl;
        cout << "gnode " << gnode << endl;
        cout << "vnode " << vnode << endl;
        cout << "gnode reference to vnode " << gnode_r << endl;
        cout << "gnode pointer to vnode " << (*gnode_p) << endl;
*/
        delete tree;
    }

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
