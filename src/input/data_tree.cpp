#include <ostream>
#include <string>
#include "../system/system.hh"
#include "data_tree.hpp"
#include "Generic_node.hpp"
#include "Record_node.hpp"
#include "Vector_node.hpp"
#include "Value_node.hpp"

#include "json_spirit.h"

#include "input_interface.hh"
#include "input_type.hh"

namespace flow {

/*!
 * @brief State machine helper for filtering comments out of input JSON stream (string).
 * @param in_char      Single input character from stream (string).
 * @param out_string   String with filtered output.
 * @param reset_state  Set reset_state to TRUE, when starting to filter new stream (string).
 */
void Data_tree::filter_stm( const char in_char, std::string & out_string, bool reset_state = false )
{
    typedef enum {GO, GO_BSL, IN_QUOTE, IN_QUOTE_BSL, IN_COMMENT, IN_COMMENT_BSL} states;
    static states state;
    states next_state;

    if ( reset_state )
    {
        state = GO;
        return;
    } else {
        next_state = state; //Failsafe for 'uninitialized variable' when programmer makes mistake.
    }

    switch (state) {
    case GO:
        switch (in_char) {
        case '#':
            next_state = IN_COMMENT;
            break;
        case '\\':
            next_state = GO_BSL;
            break;
        case '"':
            out_string.push_back(in_char);
            next_state = IN_QUOTE;
            break;
        default:
            out_string.push_back(in_char);
            next_state = GO;
            break;
        }
        break;
    case GO_BSL:
        out_string.push_back('\\');
        out_string.push_back(in_char);
        next_state = GO;
        break;
    case IN_QUOTE:
        switch (in_char) {
        case '"':
            out_string.push_back(in_char);
            next_state = GO;
            break;
        case '\\':
            next_state = IN_QUOTE_BSL;
            break;
        default:
            out_string.push_back(in_char);
            next_state = IN_QUOTE;
            break;
        }
        break;
    case IN_QUOTE_BSL:
        out_string.push_back('\\');
        out_string.push_back(in_char);
        next_state = IN_QUOTE;
        break;
    case IN_COMMENT:
        switch (in_char) {
        case '\n':
            next_state = GO;
            break;
        case '\r':
            next_state = GO;
            break;
        case '\\':
            next_state = IN_COMMENT_BSL;
            break;
        default:
            next_state = IN_COMMENT;
            break;
        }
        break;
    case IN_COMMENT_BSL:
        switch (in_char) {
        case '\n':
        case '\r':
            next_state = IN_COMMENT_BSL;
            break;
        default:
            next_state = IN_COMMENT;
            break;
        }
        break;
    default:
        xprintf(PrgErr,"Previous state is unknown, should never happen.");
        break;
    }

    state = next_state;
}


string Data_tree::flow_json_filter( std::istream& is )
{
    string ret_s;

    filter_stm( '0', ret_s, true );

    while ( is.good() )
    {
        char c;
        c = is.get();
        filter_stm( c, ret_s );
    }

    return ret_s;
}

string Data_tree::flow_json_filter( const std::string& s )
{
    string ret_s;
    int pos = 0;
    int size;

    filter_stm( '0', ret_s, true );

    size = s.size();
    while ( pos < size )
    {
        char c;
        c = s.at(pos++);
        filter_stm( c, ret_s );
    }

    return ret_s;
}

Data_tree::Data_tree( const std::string& s )
{
    //use json_spirit mValue (using map) and not Value (vector)
    //vector is exponentially slower for large data
    json_spirit::mValue json_root; //root of loaded JSON file (or stream)

    err_status = false;

    //load JSON from string
    if ( !json_spirit::read( flow_json_filter(s),json_root) )
    {
        err_status = true;
        return;
    }

    //build reference tree
    if ( !tree_build( json_root, node_head ) )
    {
        err_status = true;
        return;
    }
}

Data_tree::Data_tree( std::istream& is )
{
    //use json_spirit mValue (using map) and not Value (vector)
    //vector is exponentially slower for large data
    json_spirit::mValue json_root; //root of loaded JSON file (or stream)

    err_status = false;

    //load JSON from string
    if ( !json_spirit::read( flow_json_filter(is), json_root) )
    {
        err_status = true;
        return;
    }

    //build reference tree
    if ( !tree_build( json_root, node_head ) )
    {
        err_status = true;
        return;
    }
}

Generic_node * Data_tree::new_node( const json_spirit::mValue json_node, Generic_node * prev_node )
{
    Generic_node * gnp = NULL;

    switch (json_node.type()) {
    case json_spirit::obj_type:
        gnp = new Record_node(prev_node);
        break;
    case json_spirit::array_type:
        gnp = new Vector_node(prev_node);
        break;
    case json_spirit::str_type:
        gnp = new Value_node(prev_node, json_node.get_str());
        break;
    case json_spirit::bool_type:
        gnp = new Value_node(prev_node, json_node.get_bool());
        break;
    case json_spirit::int_type:
        gnp = new Value_node(prev_node, json_node.get_int());
        break;
    case json_spirit::real_type:
        gnp = new Value_node(prev_node, json_node.get_real());
        break;
    case json_spirit::null_type:
        gnp = new Value_node(prev_node);
        break;
    default:
        xprintf( PrgErr, "Unknown node type in original JSON tree." );
        break;
    }

    return gnp;
}

bool Data_tree::tree_build_recurse( json_spirit::mValue json_root, Generic_node & prev_node )
{
    switch ( json_root.type() ) {
    case json_spirit::obj_type: {
            json_spirit::mObject::iterator it;
            Record_node & o_node = prev_node.as_record();

            for( it = json_root.get_obj().begin(); it != json_root.get_obj().end(); ++it )
            {
                Generic_node * gnp = new_node(it->second, &prev_node);
                o_node.insert_key( it->first, gnp );

                switch (it->second.type()) {
                case json_spirit::obj_type: //Record need recursive build
                    //cout << "KEY: " << it->first << " ";
                    //cout << "Record: going deep..." << endl;
                    if (!tree_build_recurse(it->second, *gnp)) {
                        xprintf( Warn, "Recursive tree build failed at node: %s \n", it->first.c_str() );
                        return false;
                    }
                    //cout << "Record: going up..." << endl;
                    break;
                case json_spirit::array_type: //Array need recursive build
                    //cout << "KEY: " << it->first << " ";
                    //cout << "Array: going deep..." << endl;
                    if (!tree_build_recurse(it->second, *gnp)) {
                        xprintf( Warn, "Recursive tree build failed at node: %s \n", it->first.c_str() );
                        return false;
                    }
                    //cout << "Array: going up..." << endl;
                    break;
                default: //Other types of nodes do not need special treatment.
                    break;
                }
            }
        }
        break;
    case json_spirit::array_type: {
            Vector_node & v_node = prev_node.as_vector();
            for( unsigned int i = 0; i < json_root.get_array().size(); ++i )
            {
                Generic_node * gnp = new_node(json_root.get_array().at(i), & prev_node);
                v_node.insert_item( i, *gnp );

                switch (json_root.get_array().at(i).type()) {
                case json_spirit::obj_type: //Record need recursive build
                    //cout << "ID: " << i << " ";
                    //cout << "Record: going deep..." << endl;
                    if (!tree_build_recurse(json_root.get_array().at(i), *gnp)) {
                        xprintf( Warn, "Recursive tree build failed at ID: %u \n", i );
                        return false;
                    }
                    //cout << "Record: going up..." << endl;
                    break;
                case json_spirit::array_type: //Array need recursive build
                    //cout << "ID: " << i << " ";
                    //cout << "Array: going deep..." << endl;
                    if (!tree_build_recurse(json_root.get_array().at(i), *gnp)) {
                        xprintf( Warn, "Recursive tree build failed at ID: %u \n", i );
                        return false;
                    }
                    //cout << "Array: going up..." << endl;
                    break;
                default: //Other types of nodes do not need special treatment.
                    break;
                }
            }
        }
        break;
    case json_spirit::str_type: //no need to recurse for scalar value
        break;
    case json_spirit::bool_type: //no need to recurse for scalar value
        break;
    case json_spirit::int_type: //no need to recurse for scalar value
        break;
    case json_spirit::real_type: //no need to recurse for scalar value
        break;
    case json_spirit::null_type: //no need to recurse for scalar value
        break;
    default:
        xprintf( PrgErr, "Unknown type of original JSON node." );
        return false;
        break;
    }

    return true;
}

bool Data_tree::tree_build( const json_spirit::mValue json_root, Generic_node & head_node )
{
    //JSON not read OK
    if ( err_status )
        return false;

    //on highest level JSON can contain: 1 record or 1 array, nothing else
    // [...,...,...]
    // {...}

    switch (json_root.type()) {
    case json_spirit::obj_type:
        return tree_build_recurse(json_root, head_node);
        break;
    case json_spirit::array_type:
        xprintf( Warn, "Top-level element in JSON data is ARRAY - NOT IMPLEMENTED.");
        return false;
        break;
    default:
        xprintf( Warn, "Top-level element in JSON data is not RECORD or ARRAY.");
        return false;
        break;
    }

    return true;
}


bool Data_tree::refs_scandel( Generic_node & head_node, vector< tree_ref > & vrefs ) {
    return true;
}

bool Data_tree::refs_unpack( Generic_node & head_node, vector< tree_ref > & vrefs ) {
    return true;
}

bool Data_tree::refs_process(Generic_node& head_node) {

    vector< tree_ref > vrefs;

    if ( !refs_scandel( head_node, vrefs ) )
        return false;

    if ( !refs_unpack( head_node, vrefs ) )
        return false;

    return true;
}

ostream & operator<<(ostream & stream, Data_tree & tree )
{
    stream << tree.node_head;
    return stream;
}

}
