#include "../system/system.hh"
#include "Record_node.hpp"
#include "Value_node.hpp"
#include "Vector_node.hpp"

namespace flow {

Generic_node & Vector_node::get_item(const size_t id) {
    if (id >= value_array_.size()) {
        //out of range - return empty
        return *empty_node_generic_;
    } else {
        return *value_array_[id];
    }
}

Generic_node & Vector_node::get_item(const size_t id, Generic_node & default_tree)
{
    if (id >= value_array_.size()) {
        //out of range - return default
        return default_tree;
    } else {
        return *value_array_[id];
    }
}

Generic_node & Vector_node::get_item_check(const size_t id, int & err_code) {
    if (id >= value_array_.size()) {
        //out of range - return empty & error
        err_code = 1;
        return *empty_node_generic_;
    } else {
        err_code = 0;
        return *value_array_[id];
    }
}

ostream & operator <<(ostream & stream, Vector_node & node) {
    size_t i;
    size_t size;

    size = node.value_array_.size();

    stream << "[";
    for ( i = 0; i < size; ++i )
    {
        stream << node.get_item( i );
        if ( (i+1) < size )
            stream << ",";
    }
    stream << "]";
    return stream;
}
void Vector_node::insert_item_parent( const size_t id, Generic_node & node, Generic_node * parent ) {
    if (id < value_array_.size()) {
        value_array_.at(id) = &node;
    } else {
        value_array_.resize(id+1);
        value_array_.at(id) = &node;
    }
    node.set_parent_node( parent );
}

void Vector_node::insert_item( const size_t id, Generic_node & node) {
    insert_item_parent( id, node, this );
}

void Vector_node::delete_item(const size_t id) {
    if (id < value_array_.size()) {
        delete_node_unpacked( id );
        value_array_.erase( value_array_.begin() + id );
    }
}

void Vector_node::delete_node_unpacked( const size_t id ) {
    value_array_[id]->~Generic_node();
    delete value_array_[id];
    value_array_[id] = NULL;
}

Vector_node::Vector_node( const Vector_node & to_copy ) {
    size_t size = to_copy.value_array_.size();
    value_array_.clear();
    value_array_.resize(size);
    size_t i;

    for ( i = 0; i < size; ++i )
    {
        switch ( to_copy.value_array_[i]->get_type() ) {
        case type_string:
        case type_number:
        case type_bool:
        case type_null:
            insert_item( i, * new Value_node( to_copy.value_array_[i]->as_value() ) );
            break;
        case type_record:
            insert_item( i, * new Record_node( to_copy.value_array_[i]->as_record() ) );
            break;
        case type_vector:
            insert_item( i, * new Vector_node( to_copy.value_array_[i]->as_vector() ) );
            break;
        default:
            xprintf( PrgErr, "Err: instance of Generic_node should not exist." );
            break;
        }
    }
}

Vector_node::~Vector_node() {
    size_t size = value_array_.size();
    size_t i;

    for ( i = 0; i < size; ++i )
    {
        delete_node_unpacked( i );
    }
    value_array_.clear();
}

} //namespace
