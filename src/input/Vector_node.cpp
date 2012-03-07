#include "Vector_node.hpp"

namespace flow {

Generic_node & Vector_node::get_item(const size_t id) {
    if (id >= value_array_.size()) {
        //mimo rozsah pole, vratime prazdnou
        return *empty_node_generic_;
    } else {
        return *value_array_[id];
    }
}

Generic_node & Vector_node::get_item(const size_t id, Generic_node & default_tree)
{
    if (id >= value_array_.size()) {
        //mimo rozsah pole, vratime default
        return default_tree;
    } else {
        return *value_array_[id];
    }
}

Generic_node & Vector_node::get_item_check(const size_t id, int & err_code) {
    if (id >= value_array_.size()) {
        //mimo rozsah pole, vratime empty, s chybou
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
        cout << node.get_item( i );
        if ( (i+1) < size )
            cout << ",";
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
        {
            //TODO: nejakou lepsi hlasku...
            cout<<"Array FUCK"<<endl;
        }

        value_array_.push_back( &node );
    }
}

Vector_node::~Vector_node() {
    //TODO: deep destructor?
}

}
