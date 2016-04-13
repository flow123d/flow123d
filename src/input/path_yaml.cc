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
 * @file    path_yaml.cc
 * @brief   
 */

#include "input/path_yaml.hh"
#include "input/reader_to_storage.hh"

#include "system/global_defs.h"
#include "system/system.hh"


namespace Input {
using namespace std;


PathYAML::PathYAML(istream &in)
: PathBase()
{
    nodes_.push_back( YAML::Load( in ) );

    json_type_names.push_back("YAML map");
    json_type_names.push_back("YAML sequence");
    json_type_names.push_back("YAML string");
    json_type_names.push_back("YAML bool");
    json_type_names.push_back("YAML int");
    json_type_names.push_back("YAML real");
    json_type_names.push_back("YAML null");
    json_type_names.push_back("other scalar type");
    json_type_names.push_back("undefined type");
}

PathYAML::~PathYAML()
{
    this->go_to_root();
}


bool PathYAML::down(unsigned int index) {
	FEAL_DEBUG_ASSERT(head().IsSequence()).error("Head node must be of type Array.");

    if ( index >= head().size() ) return false;
    path_.push_back( make_pair( index, string("") ) );
    nodes_.push_back( head()[index] );

    return true;
}


bool PathYAML::down(const string& key) {
	FEAL_DEBUG_ASSERT(head().IsMap()).error("Head node must be of type Record.");

    if ( head()[key] ) {
    	path_.push_back( make_pair( (int)(-1), key) );
    	nodes_.push_back( head()[key] );
    } else {
        return false;
    }
    return true;
}


void PathYAML::up() {
    if (path_.size() > 1) {
        path_.pop_back();
        nodes_.pop_back();
    }
}


bool PathYAML::is_null_type() const {
	return head().IsNull();
}


bool PathYAML::get_bool_value() const {
	if (head().IsScalar()) {
		try {
			return head().as<bool>();
		} catch (YAML::Exception) {
	        THROW( ReaderToStorage::ExcInputError() );
		}
	} else {
        THROW( ReaderToStorage::ExcInputError() );
	}
	return false;
}


std::int64_t PathYAML::get_int_value() const {
	if (head().IsScalar()) {
		try {
			return head().as<std::int64_t>();
		} catch (YAML::Exception) {
	        THROW( ReaderToStorage::ExcInputError() );
		}
	} else {
        THROW( ReaderToStorage::ExcInputError() );
	}
	return 0;
}


double PathYAML::get_double_value() const {
	if (head().IsScalar()) {
		try {
			return head().as<double>();
		} catch (YAML::Exception) {
	        THROW( ReaderToStorage::ExcInputError() );
		}
	} else {
        THROW( ReaderToStorage::ExcInputError() );
	}
	return 0.0;
}


std::string PathYAML::get_string_value() const {
	if (head().IsScalar()) {
		try {
			return head().as<std::string>();
		} catch (YAML::Exception) {
	        THROW( ReaderToStorage::ExcInputError() );
		}
	} else {
        THROW( ReaderToStorage::ExcInputError() );
	}
	return "";
}


unsigned int PathYAML::get_node_type_index() const {
	switch (head().Type()) {
	  case YAML::NodeType::Null: return ValueTypes::null_type;
	  case YAML::NodeType::Scalar: return ValueTypes::scalar_type;
	  case YAML::NodeType::Sequence: return ValueTypes::array_type;
	  case YAML::NodeType::Map: return ValueTypes::obj_type;
	  default: return ValueTypes::undef_type;
	}
}


bool PathYAML::get_record_key_set(std::set<std::string> &keys_list) const {
	if ( head().IsMap() ) {
		for (YAML::const_iterator it=head().begin(); it!=head().end(); ++it) {
			keys_list.insert( it->first.as<std::string>() );
		}
        return true;
    }

	return false;
}


int PathYAML::get_array_size() const {
	if (head().IsSequence()) {
		return head().size();
	}
	return -1;
}


bool PathYAML::is_record_type() const {
	return head().IsMap();
}


bool PathYAML::is_array_type() const {
	return head().IsSequence();
}


PathYAML * PathYAML::clone() const {
	return new PathYAML(*this);
}


PathBase * PathYAML::find_ref_node()
{
    return NULL;
}



std::string PathYAML::get_descendant_name() const {
	std::string tag = head().Tag();
	if (tag == "?") {
		return "";
	} else {
		return tag.erase(0, 1); // tag starts with '!' char
	}
}



std::ostream& operator<<(std::ostream& stream, const PathYAML& path) {
    path.output(stream);
    return stream;
}


} // namespace Input
