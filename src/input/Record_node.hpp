#ifndef RECORD_NODE_HPP_
#define RECORD_NODE_HPP_

#include <string>
#include "Generic_node.hpp"

namespace flow {

/*!
 * @brief Record node - represents JSON construct "{}"
 *        with pairs of "string" : any_node
 */
class Record_node: public Generic_node {
    map< string, Generic_node & > record_;
public:
    Record_node():Generic_node(type_record) {}
    Record_node( Generic_node & prev_node ):Generic_node(type_record, prev_node) {}

    virtual Generic_node & get_item( const int id ) {
        //Vector-like access to Record - return empty
        return *empty_node_generic_;
    }
    virtual Generic_node & get_item( const size_t id, Generic_node & default_tree ) {
        //Vector-like access to Record - return empty
        return default_tree;
    }
    virtual Generic_node & get_item_check( const size_t id, int & err_code ) {
        //Vector-like access to Record - return empty & error
        err_code = 1;
        return *empty_node_generic_;
    }

    //insert new pair of "key":any_node into record
    void insert_key( const string & key, Generic_node & node );

    virtual Generic_node & get_key( const string & key );
    virtual Generic_node & get_key( const string & key, Generic_node & default_tree );
    virtual Generic_node & get_key_check( const string & key, int & err_code );

    virtual Record_node & as_record( void ) { return (*this); }

    friend ostream & operator<<( ostream & stream, Record_node & node );

    virtual ~Record_node();
};

}

#endif /* RECORD_NODE_HPP_ */
