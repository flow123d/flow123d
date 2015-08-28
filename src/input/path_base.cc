/*
 * path_base.cc
 *
 *  Created on: May 7, 2012
 *      Author: jb
 */


#include "input/path_base.hh"
#include "input/input_exception.hh"


namespace Input {
using namespace std;


PathBase::PathBase() {
	path_.push_back( make_pair( (int)(-1), string("/") ) );
}


void PathBase::output(ostream &stream) const {
    if (level() == 0) {
        stream << "/";
        return;
    }

    for(vector< pair<int, string> >::const_iterator it = path_.begin()+1; it != path_.end(); ++it) {
        if ( it->first < 0 ) {
            stream << "/" << it->second;
        } else {
            stream << "/" << it->first;
        }
    }
}



string PathBase::as_string() const {
    stringstream ss;
    output(ss);
    return ss.str();
}



void PathBase::go_to_root() {
	while (path_.size() > 1) {
		this->up();
	}
}



std::string PathBase::get_node_type(unsigned int type_idx) const {
	return json_type_names[ type_idx ];
}



} // namespace Input
