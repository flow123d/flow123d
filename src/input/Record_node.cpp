#include "Record_node.hpp"

namespace flow {

Generic_node & Record_node::get_key(const string & key) {
    map<string, Generic_node &>::iterator it;

    it = record_.find(key);
    if (it == record_.end()) {
        return *empty_node_generic_;
    } else {
        return it->second;
    }
}

Generic_node & Record_node::get_key(const string & key, Generic_node & default_tree) {
    map<string, Generic_node &>::iterator it;

    it = record_.find(key);
    if (it == record_.end()) {
        return default_tree;
    } else {
        return it->second;
    }
}

ostream & operator <<(ostream & stream, Record_node & node) {
    map< string, Generic_node & >::iterator it;

    stream << "{";

    it=node.record_.begin();
    while ( it != node.record_.end() )
    {
        stream << "\"" << it->first << "\":" << it->second;
        ++it;
        if ( it != node.record_.end() )
            stream << ",";
    }

    stream << "}";
    return stream;
}

Generic_node & Record_node::get_key_check(const string & key, int & err_code) {
    map<string, Generic_node &>::iterator it;

    it = record_.find(key);
    if (it == record_.end()) {
        err_code = 1;
        return *empty_node_generic_;
    } else {
        err_code = 0;
        return it->second;
    }
}

void Record_node::delete_key(const string& key) {
    map<string, Generic_node &>::iterator it;

    it = record_.find(key);
    if (it != record_.end()) {
        record_.erase(it);
    }
}

void Record_node::insert_key_parent( const string & key, Generic_node & node, Generic_node * parent ) {
    delete_key( key );
    record_.insert( pair<string,Generic_node&>(key,node) );
    node.set_parent_node( parent );
}

void Record_node::insert_key(const string & key, Generic_node & node) {
    insert_key_parent( key, node, this );
}

Record_node::~Record_node() {
    //TODO deep destructor

}

}
