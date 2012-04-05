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

void Vector_node::insert_item( const size_t id, Generic_node & node) {
    if (id < value_array_.size()) {
        value_array_.at(id) = &node;
    } else {
        //have room?
        if ( value_array_.size() == value_array_.max_size())
            xprintf( PrgErr, "Memory allocation error." );

        value_array_.push_back( &node );
    }
}

Vector_node::~Vector_node() {
    //call proper destructor for all types
    size_t size = value_array_.size();
    size_t i;

    for ( i = 0; i < size; ++i )
    {
        switch ( value_array_[i]->get_type() ) {
        case type_record:
            dynamic_cast < Record_node * > (value_array_[i])->~Record_node();
            break;
        case type_vector:
            dynamic_cast < Vector_node * > (value_array_[i])->~Vector_node();
            break;
        case type_null:
        case type_bool:
        case type_number:
        case type_string:
            dynamic_cast < Value_node * > (value_array_[i])->~Value_node();
            break;
        case type_generic:
        default:
        break;
        }
        delete value_array_[i];
        value_array_[i] = NULL;
    }
    value_array_.clear();
}

} //namespace
