/*!
 *
﻿ * Copyright (C) 2015 Technical University of Liberec.  All rights reserved.
 * 
 * This program is free software; you can redistribute it and/or modify it under
 * the terms of the GNU General Public License version 3 as published by the
 * Free Software Foundation. (http://www.gnu.org/licenses/gpl-3.0.en.html)
 * 
 * This program is distributed in the hope that it will be useful, but WITHOUT
 * ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 * FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.
 *
 * 
 * @file    path_base.cc
 * @brief   
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
