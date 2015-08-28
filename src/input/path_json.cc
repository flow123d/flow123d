/*
 * path_json.cc
 *
 *  Created on: May 7, 2012
 *      Author: jb
 */



#include "input/path_json.hh"
#include "input/reader_to_storage.hh"
#include "input/comment_filter.hh"

#include "json_spirit/json_spirit_error_position.h"
#include <boost/iostreams/device/file.hpp>
#include <boost/iostreams/filtering_stream.hpp>


namespace Input {
using namespace std;


PathJSON::PathJSON(istream &in)
: PathJSON()
{
    io::filtering_istream filter_in;

    filter_in.push(uncommenting_filter());
    filter_in.push(in);

    Node root_node;

    try {
        json_spirit::read_or_throw( filter_in, root_node);
    } catch (json_spirit::Error_position &e ) {
        THROW( ReaderToStorage::ExcNotJSONFormat() << ReaderToStorage::EI_JSONLine(e.line_) << ReaderToStorage::EI_JSONColumn(e.column_)
        	<< ReaderToStorage::EI_JSONReason(e.reason_));
    }

    nodes_.push_back( new Node(root_node) );
}


PathJSON::PathJSON()
: PathBase()
{
    /* from json_spirit_value.hh:
     * enum Value_type{ obj_type, array_type, str_type, bool_type, int_type, real_type, null_type };
     */
    json_type_names.push_back("JSON object");
    json_type_names.push_back("JSON array");
    json_type_names.push_back("JSON string");
    json_type_names.push_back("JSON bool");
    json_type_names.push_back("JSON int");
    json_type_names.push_back("JSON real");
    json_type_names.push_back("JSON null");
    json_type_names.push_back(""); //scalar type
    json_type_names.push_back(""); //undefined type
}


bool PathJSON::down(unsigned int index)
{
    const Node * head_node = nodes_.back();
    const json_spirit::mArray &array = head_node->get_array(); // the type should be checked in make_storage

    if (  index >= array.size()) return false;
    path_.push_back( make_pair( index, string("") ) );
    nodes_.push_back( &( array[index] ) );
    return true;
}


bool PathJSON::down(const string& key)
{
    const Node * head_node = nodes_.back();
    const json_spirit::mObject &obj = head_node->get_obj(); // the type should be checked in make_storage

    json_spirit::mObject::const_iterator it = obj.find(key);
    if (it == obj.end()) {
        return false;
    } else {
        path_.push_back( make_pair( (int)(-1), key) );
        nodes_.push_back( &( it->second ) );
    }
    return true;
}


void PathJSON::up()
{
    if (path_.size() > 1) {
        path_.pop_back();
        nodes_.pop_back();
    }
}


/**
 * This returns path to reference given by address in ref_address.
 *
 */
PathBase * PathJSON::find_ref_node()
{
    namespace ba = boost::algorithm;

    const Node * head_node = nodes_.back();
    if (head_node->type() != json_spirit::obj_type) return NULL;
    const json_spirit::mObject &obj = head_node->get_obj();
    if (obj.size() != 1) return NULL;
    if (obj.begin()->first != "REF") return NULL;

    const Node &ref_node = obj.begin()->second;
    if (ref_node.type() != json_spirit::str_type) {
        THROW( ExcRefOfWrongType() << EI_ErrorAddress(this->as_string()) );

    }
    string ref_address = ref_node.get_str();

    PathJSON * ref_path = this->clone();
    string::size_type pos = 0;
    string::size_type new_pos = 0;
    string address = ref_address + '/';
    string tmp_str;
    bool relative_ref = false;

    std::set<string>::iterator it = previous_references_.find(ref_address);
    if (it == previous_references_.end()) {
    	ref_path->previous_references_.insert(ref_address);
    } else {
    	THROW( ExcReferenceNotFound() << EI_RefAddress(this->as_string()) << EI_ErrorAddress(ref_path->as_string()) << EI_RefStr(ref_address)
    	       << EI_Specification("cannot follow reference") );
    }

    while ( ( new_pos=address.find('/',pos) ) != string::npos ) {
        tmp_str = address.substr(pos, new_pos - pos);
        if (pos==0 && tmp_str == "") {
            // absolute path
            ref_path->go_to_root();

        } else if ( ba::all( tmp_str, ba::is_digit()) ) {
            // integer == index in array
            if ( !ref_path->is_array_type() ) {
                THROW( ExcReferenceNotFound() << EI_RefAddress(this->as_string()) << EI_ErrorAddress(ref_path->as_string()) << EI_RefStr(ref_address)
                        << EI_Specification("there should be Array") );
            }

            if ( !ref_path->down( atoi(tmp_str.c_str()) ) ) {
                THROW( ExcReferenceNotFound() << EI_RefAddress(this->as_string()) << EI_ErrorAddress(ref_path->as_string()) << EI_RefStr(ref_address)
                        << EI_Specification("index out of size of Array") );
            }

        } else if (tmp_str == "..") {
        	relative_ref = true;
            if (ref_path->level() <= 0 ) {
                THROW( ExcReferenceNotFound() << EI_RefAddress(this->as_string()) << EI_ErrorAddress(ref_path->as_string()) << EI_RefStr(ref_address)
                        << EI_Specification("can not go up from root") );
            }
            ref_path->up();

        } else {
            if ( !ref_path->is_record_type() )
                THROW( ExcReferenceNotFound() << EI_RefAddress(this->as_string()) << EI_ErrorAddress(ref_path->as_string()) << EI_RefStr(ref_address)
                        << EI_Specification("there should be Record") );
            if ( !ref_path->down(tmp_str) )
                THROW( ExcReferenceNotFound() << EI_RefAddress(this->as_string()) << EI_ErrorAddress(ref_path->as_string()) << EI_RefStr(ref_address)
                        << EI_Specification("key '"+tmp_str+"' not found") );
        }
        pos = new_pos+1;
    }
    if (relative_ref) {
    	xprintf(Msg, "Referencing '%s' to '%s'.\n", this->as_string().c_str(), ref_path->as_string().c_str());
    }
    return ref_path;
}



bool PathJSON::is_null_type() const {
	return head()->type() == json_spirit::null_type;
}


bool PathJSON::get_bool_value() const {
    if (head()->type() == json_spirit::bool_type) {
        return head()->get_bool();
    } else {
        THROW( ReaderToStorage::ExcInputError()  );
    }
	return false;
}



std::int64_t PathJSON::get_int_value() const {
    if (head()->type() == json_spirit::int_type) {
        return head()->get_int64();
    } else {
        THROW( ReaderToStorage::ExcInputError() );
    }
	return 0;
}



double PathJSON::get_double_value() const {
    auto value_type = head()->type();
    if (    value_type== json_spirit::real_type
         || value_type == json_spirit::int_type) {
        return head()->get_real();
    } else {
        THROW( ReaderToStorage::ExcInputError() );
    }
	return 0.0;
}



std::string PathJSON::get_string_value() const {
    if (head()->type() == json_spirit::str_type) {
        return head()->get_str();
    } else {
        THROW( ReaderToStorage::ExcInputError() );
    }
	return "";
}



unsigned int PathJSON::get_node_type_index() const {
	return head()->type();
}



bool PathJSON::get_record_key_set(std::set<std::string> &keys_list) const {
	if ( this->is_record_type() ) {
    	const json_spirit::mObject & j_map = head()->get_obj();
        json_spirit::mObject::const_iterator map_it;
        std::set<string>::iterator set_it;

        for( map_it = j_map.begin(); map_it != j_map.end(); ++map_it) {
        	keys_list.insert(map_it->first);
        }

        return true;
    }

	return false;
}



int PathJSON::get_array_size() const {
    if (head()->type() == json_spirit::array_type) {
        const json_spirit::mArray & j_array = head()->get_array();
        return j_array.size();
    }

    return -1;
}



bool PathJSON::is_record_type() const {
	return head()->type() == json_spirit::obj_type;
}



bool PathJSON::is_array_type() const {
	return head()->type() == json_spirit::array_type;
}



PathJSON * PathJSON::clone() const {
	return new PathJSON(*this);
}



std::string PathJSON::get_descendant_name() const {
	std::string desc_name = "";
	PathJSON type_path(*this);
	if ( type_path.down("TYPE") ) {
	    //check if TYPE key is reference
		PathBase * ref_path = type_path.find_ref_node();
	    if (ref_path) {
	        desc_name = ref_path->get_string_value();
	    	delete ref_path;
	    } else {
	    	desc_name = type_path.get_string_value();
		}
	}

	return desc_name;
}



void PathJSON::remember_reference() {
	previous_references_.insert(as_string());
}



std::ostream& operator<<(std::ostream& stream, const PathJSON& path) {
    path.output(stream);
    return stream;
}


} // namespace Input
