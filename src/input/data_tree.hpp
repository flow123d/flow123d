#ifndef DATA_TREE_HPP_
#define DATA_TREE_HPP_

#include <iostream>
#include <string>

#include "Generic_node.hpp"
#include "Record_node.hpp"
#include "Vector_node.hpp"
#include "Value_node.hpp"

// TODO: try to remove this from header file, should be possible by forwar declaration
#include "json_spirit.h"


namespace flow {

/*!
 * @brief Manages loading JSON file and building node tree,
 *        that contains loaded data. Provides various access methods to loaded data.
 *
 */
class Data_tree {
protected:
    Record_node node_head;         //root of node hierarchy
private:
    //TODO: Copy constructor? (Deep copy?)
    Data_tree():err_status(false) {};

    Generic_node * new_node( const json_spirit::mValue json_node, Generic_node & prev_node );
    bool tree_build( const json_spirit::mValue json_root, Generic_node & head_node );
    bool tree_build_recurse( json_spirit::mValue json_root, Generic_node & node );
    void filter_stm( const char in_char, string & out_string, bool reset_state );
public:
    bool err_status;

    string flow_json_filter( const std::string& s ); //< input filter for string
    string flow_json_filter( std::istream& is ); //< input filter for stream

    Data_tree( const std::string& s ); //< build tree from JSON in string
    Data_tree( std::istream& is );     //< build tree from JSON in stream

    friend ostream & operator<<( ostream & stream, Data_tree & tree );

    Generic_node & get_head() { return node_head; }
};


}

#endif /* DATA_TREE_HPP_ */
