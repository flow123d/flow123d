#include "../system/system.hh"
#include "Record_node.hpp"
#include "Vector_node.hpp"
#include "Value_node.hpp"

namespace flow {

Generic_node & Record_node::get_key(const string & key) {
    map< string, Generic_node * >::iterator it;

    it = record_.find(key);
    if (it == record_.end()) {
        return *empty_node_generic_;
    } else {
        return *(it->second);
    }
}

Generic_node & Record_node::get_key(const string & key, Generic_node & default_tree) {
    map< string, Generic_node * >::iterator it;

    it = record_.find(key);
    if (it == record_.end()) {
        return default_tree;
    } else {
        return *(it->second);
    }
}

Generic_node & Record_node::get_key_check(const string & key, int & err_code) {
    map< string, Generic_node * >::iterator it;

    it = record_.find(key);
    if (it == record_.end()) {
        err_code = 1;
        return *empty_node_generic_;
    } else {
        err_code = 0;
        return *(it->second);
    }
}

void Record_node::delete_key(const string& key) {
    map< string, Generic_node * >::iterator it;

    it = record_.find(key);
    if (it != record_.end()) {
        it->second->~Generic_node();
        record_.erase(it);
    }
}

void Record_node::insert_key_parent( const string & key, Generic_node * node, Generic_node * parent ) {
    delete_key( key );
    record_.insert( pair< string, Generic_node * >(key, node) );
    node->set_parent_node( parent );
}

void Record_node::insert_key(const string & key, Generic_node * node) {
    insert_key_parent( key, node, this );
}

ostream & operator <<(ostream & stream, Record_node & node) {
    map< string, Generic_node * >::iterator it;

    stream << "{";

    it=node.record_.begin();
    while ( it != node.record_.end() )
    {
        stream << "\"" << it->first << "\":" << *(it->second);
        ++it;
        if ( it != node.record_.end() )
            stream << ",";
    }

    stream << "}";
    return stream;
}

Record_node::Record_node( Record_node & to_copy ) {
    map< string, Generic_node * >::iterator it;

    it = to_copy.record_.begin();
    while ( it != to_copy.record_.end() ) {
        switch ( it->second->get_type() ) {
        case type_string:
        case type_number:
        case type_bool:
        case type_null:
            insert_key( it->first, new Value_node( it->second->as_value() ));
            break;
        case type_record:
            insert_key( it->first, new Record_node( it->second->as_record() ));
            break;
        case type_vector:
            insert_key( it->first, new Vector_node( it->second->as_vector() ));
            break;
        default:
            xprintf( PrgErr, "Err: instance of Generic_node should not exist." );
            break;
        }
        ++it;
    }
}

Record_node::~Record_node() {
    map< string, Generic_node * >::iterator it;
    while ( record_.size() > 0 ) {
        it = record_.begin();
        it->second->~Generic_node();
        record_.erase(it);
    }
}

} //namespace
