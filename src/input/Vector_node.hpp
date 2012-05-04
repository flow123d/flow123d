#ifndef VECTOR_NODE_HPP_
#define VECTOR_NODE_HPP_

#include "Generic_node.hpp"

namespace flow {

/*!
 * @brief Array node - represents JSON construct "[]",
 *        that is heterogeneous vector of any JSON construct
 */
class Vector_node : public Generic_node {
    //Vector of references not allowed - need to use pointers.
    vector< Generic_node * >  value_array_;

    void delete_node_unpacked( const size_t id );

public:
    Vector_node ():Generic_node(type_vector)                                      {}
    Vector_node ( Generic_node * prev_node ):Generic_node(type_vector, prev_node) {}
    Vector_node ( const Vector_node & to_copy );

    virtual Generic_node & get_item( const size_t id );
    virtual Generic_node & get_item( const size_t id, Generic_node & default_tree );
    virtual Generic_node & get_item_check( const size_t id, int & err_code );
    virtual Generic_node & get_key( const string & key ) {
        //Record-like access to Vector => return empty
        return *empty_node_generic_;
    }
    virtual Generic_node & get_key( const string & key, Generic_node & default_tree ) {
        //Record-like access to Vector => return default
        return default_tree;
    }
    virtual Generic_node & get_key_check( const string & key, int & err_code ) {
        //Record-like access to Vector => return empty & error
        err_code = 1;
        return *empty_node_generic_;
    }

    void insert_item( const size_t id, Generic_node & node );
    void insert_item_parent( const size_t id, Generic_node & node, Generic_node * parent );
    virtual void delete_item(const size_t id);
    virtual size_t get_array_size( void ) { return value_array_.size(); }

    virtual Vector_node & as_vector( void ) { return (*this); }
    friend ostream & operator<<( ostream & stream, Vector_node & node );

    virtual ~Vector_node();
};

}

#endif /* VECTOR_NODE_HPP_ */
