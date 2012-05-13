/**
 * Algoritmus:
 * - vse na jeden pruchod Json stromem
 * - json_to_storage( json_spirit::mvalue  json_node, Type::TypeBase tb, Storage::Node out_node)
 *   1) if-else  like swith podle typu tb
 *      volani patricne funkce
 *
 * json_spirit_reader - research:
 *  1) pouzit klidne read_or_throw variantu, beztak vola read_or_throw uvnitr
 *  2) XX = Value, Object, Array - implementation for ASCII input, using vector of pairs for Object
 *     mXX  - ASCII, use map for Object
 *     wXX  - Unicode, use vector for Object
 *     wmXX - Unicode, use map
 *  3) internally it use boost::variant for nodes of the tree
 *
 *
 *
 *
 * TODO:
 *  - Input::Storage interface to json_spirit::mvalue ...
 *    - test input
 *    -
 *  - pass through json tree check correctness, transform Object to Array
 *  - include expansion of references
 *
 *  - use stringstream instead of separate flow_json_filter( const std::string& s );
 *  - make flow_json_filter true stream filter, do not process whole string
 */

#ifndef DATA_TREE_HPP_
#define DATA_TREE_HPP_

#include <iostream>
#include <string>

#include "Generic_node.hpp"
#include "Record_node.hpp"
#include "Vector_node.hpp"
#include "Value_node.hpp"

#include "json_spirit.h" //forward declaration? too many templates for me :-/

namespace flow {

/*!
 * @brief Manages loading JSON file and building node tree,
 *        that contains loaded data. Provides various access methods to loaded data.
 *
 */
struct tree_ref {
    Generic_node * src;
    Generic_node * dst;
};

class Data_tree {
protected:
    Record_node node_head;         //root of node hierarchy
private:
    //TODO: Copy constructor? (Deep copy?)
    Data_tree():err_status(false) {};

    /**
     * For given JSON node @p json_node and constructs and returns Generic_node with appropriate type.
     *
     * It needs prev_node for back refferences, but possibly can be removed.
     */
    Generic_node * new_node( const json_spirit::mValue json_node, Generic_node * prev_node );

    bool tree_build( const json_spirit::mValue json_root, Generic_node & head_node );
    bool tree_build_recurse( json_spirit::mValue json_root, Generic_node & node );

                     /**
                      *  helper method (implements state machine for flow_json_filter methods
                      */
    void filter_stm( const char in_char, string & out_string, bool reset_state );

    /**
     * References not implemented yet !!!
     */
    bool refs_process( Generic_node & head_node );
    bool refs_scandel( Generic_node & head_node, vector< tree_ref > & vrefs );
    bool refs_unpack( Generic_node & head_node, vector< tree_ref > & vrefs );

public:
    bool err_status;

    /**
     * These methods reads stream or string and produce string with filtered comments.
     */
    string flow_json_filter( const std::string& s ); //< input filter for string
    string flow_json_filter( std::istream& is ); //< input filter for stream

    Data_tree( const std::string& s ); //< build tree from JSON in string
    Data_tree( std::istream& is );     //< build tree from JSON in stream

    friend ostream & operator<<( ostream & stream, Data_tree & tree );

    Generic_node & get_head() { return node_head; }
};


/**
 * This class holds pointer to an json_spirit::mvalue and also pointers to all its parent nodes. This
 * it can return
 */
#ifdef NOTHING
class JSONPath {
                        /**
                         * Checks if current node is a reference. If not returns and empty JSINIter, else
                         * try to find json_spirit node given by address stored in a reference.
                         */
    JSONPath refered_node();
                        /**
                         * Returns true for an empty iterator.
                         */

    bool is_null();

    deeper(const json_spirit::mvalue)
private:
    vector<const json_spirit::mvalue> path;
};
#endif

} // namespace flow


#endif /* DATA_TREE_HPP_ */
