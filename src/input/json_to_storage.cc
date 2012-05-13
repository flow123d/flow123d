/*
 * json_to_storage.cc
 *
 *  Created on: May 7, 2012
 *      Author: jb
 */

#include <boost/iostreams/device/file.hpp>
#include <boost/iostreams/filtering_stream.hpp>

#include "json_to_storage.hh"
#include "input/comment_filter.hh"
namespace Input {
using namespace std;
using namespace internal;

/********************************************
 * Implementation of public part of JSONToStorage
 */

JSONToStorage::JSONToStorage()
:storage_(NULL), root_type_(NULL)
{}


void JSONToStorage::read_stream(istream &in, const Type::TypeBase &root_type) {
    namespace io = boost::iostreams;

    io::filtering_istream filter_in;

    filter_in.push(uncommenting_filter());
    filter_in.push(in);

    JSONPath::Node node;
    json_spirit::read_or_throw( filter_in, node);

    JSONPath root_path(node);

    storage_ = make_storage(root_path, root_type_);
    root_type_ = &root_type;
}


/********************************************
 * Implementation of  internal::JSONPath
 */


JSONPath::JSONPath(const Node& root_node)
{
    path_.push_back( make_pair( (int)(-1), string("/") ) );
    nodes_.push_back( &root_node );
}

const JSONPath::Node * JSONPath::down(unsigned int index)
{
    const Node * head_node = nodes_.back();
    const json_spirit::mArray &array = head_node->get_array(); // the type should be checked in make_storage

    if (  index >= array.size()) return NULL;
    path_.push_back( make_pair( index, string("") ) );
    nodes_.push_back( &( array[index] ) );

    return nodes_.back();
}


const JSONPath::Node * JSONPath::down(const string& key)
{
    const Node * head_node = nodes_.back();
    const json_spirit::mObject &obj = head_node->get_obj(); // the type should be checked in make_storage

    json_spirit::mObject::const_iterator it = obj.find(key);
    if (it == obj.end()) {
        return NULL;
    } else {
        path_.push_back( make_pair( (int)(-1), key) );
        nodes_.push_back( &( it->second ) );
    }
    return nodes_.back();
}


void JSONPath::up()
{
    if (path_.size() > 1) {
        path_.pop_back();
        nodes_.pop_back();
    }
}

void JSONPath::go_to_root() {
    path_.resize(1);
    nodes_.resize(1);
}


bool JSONPath::get_ref_from_head(string & ref_address)
{
    const Node * head_node = nodes_.back();
    if (head_node->type() != json_spirit::obj_type) return false;
    const json_spirit::mObject &obj = head_node->get_obj();
    if (obj.size() != 1) return false;
    if (obj.begin()->first != "REF") return false;

    const Node &ref_node = obj.begin()->second;
    if (ref_node.type() != json_spirit::str_type) {
        THROW( ExcRefOfWrongType() << EI_ErrorAddress(*this) );

    }
    ref_address = ref_node.get_str();
    return true;
}

/**
 * This moves path to given reference.
 *
 * TODO:
 * << operator
 * better error messages, we need both starting path and current path
 *
 * use boost::spirit to parse address and make semantic actions.
 */
JSONPath JSONPath::find_ref_node(const string& ref_address)
{
    namespace ba = boost::algorithm;

    JSONPath ref_path(*this);

    string::size_type pos = 0;
    string::size_type new_pos = 0;
    string address = ref_address + '/';
    string tmp_str;

    while ( ( new_pos=address.find('/',pos) ) != string::npos ) {
        tmp_str = address.substr(pos, new_pos - pos);
        // DBGMSG("adr: '%s' tstr '%s' pos:%d npos:%d\n", address.c_str(), tmp_str.c_str(), pos, new_pos  );
        if (pos==0 && tmp_str == "") {
            // absolute path
            ref_path.go_to_root();

        } else if ( ba::all( tmp_str, ba::is_digit()) ) {
            // integer == index in array
            if (ref_path.head() -> type() != json_spirit::array_type) {
                THROW( ExcReferenceNotFound() << EI_RefAddress(*this) << EI_ErrorAddress(ref_path) << EI_RefStr(ref_address)
                        << EI_Specification("there should be Array") );
            }

            if ( ref_path.down( atoi(tmp_str.c_str()) ) == NULL) {
                THROW( ExcReferenceNotFound() << EI_RefAddress(*this) << EI_ErrorAddress(ref_path) << EI_RefStr(ref_address)
                        << EI_Specification("index out of size of Array") );
            }

        } else if (tmp_str == "..") {
            if (ref_path.level() <= 0 ) {
                THROW( ExcReferenceNotFound() << EI_RefAddress(*this) << EI_ErrorAddress(ref_path) << EI_RefStr(ref_address)
                        << EI_Specification("can not go up from root") );
            }
            ref_path.up();

        } else {
            if (ref_path.head() -> type() != json_spirit::obj_type)
                THROW( ExcReferenceNotFound() << EI_RefAddress(*this) << EI_ErrorAddress(ref_path) << EI_RefStr(ref_address)
                        << EI_Specification("there should be Record") );
            if ( ref_path.down(tmp_str) == NULL )
                THROW( ExcReferenceNotFound() << EI_RefAddress(*this) << EI_ErrorAddress(ref_path) << EI_RefStr(ref_address)
                        << EI_Specification("key '"+tmp_str+"' not found") );
        }
        pos = new_pos+1;
    }
    return ref_path;
}

void JSONPath::output(ostream &stream) const {
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



string JSONPath::str() {
    stringstream ss;
    output(ss);
    return ss.str();
}



std::ostream& operator<<(std::ostream& stream, const JSONPath& path) {
    path.output(stream);
    return stream;
}


/********************************************
 * Implementation of private part of JSONToStorage - make_storage dispatch
 */

StorageBase * JSONToStorage::make_storage(JSONPath &p, const Type::TypeBase *type)
{
    // first check reference
    string ref_address;
    if (p.get_ref_from_head(ref_address)) {
        // todo: mark passed references and check cyclic references

        // dereference and take data from there
        JSONPath ref_path = p.find_ref_node(ref_address);
        return make_storage( ref_path, type );
    }

    // dispatch types
    if (typeid(*type) == typeid(Type::Record)) {
        return make_storage(p, static_cast<const Type::Record *>(type) );
    } else
    if (typeid(*type) == typeid(Type::AbstractRecord)) {
        return make_storage(p, static_cast<const Type::AbstractRecord *>(type) );
    } else
    if (typeid(*type) == typeid(Type::Array)) {
        return make_storage(p, static_cast<const Type::Array *>(type) );
    } else
    if (typeid(*type) == typeid(Type::Integer)) {
        return make_storage(p, static_cast<const Type::Integer *>(type) );
    } else
    if (typeid(*type) == typeid(Type::Double)) {
        return make_storage(p, static_cast<const Type::Double *>(type) );
    } else
    if (typeid(*type) == typeid(Type::Bool)) {
        return make_storage(p, static_cast<const Type::Bool *>(type) );
    } else {
        const Type::String * string_type = dynamic_cast<const Type::String *>(type);
        if (string_type != NULL ) return make_storage(p, string_type );

        const Type::SelectionBase * selection_type = dynamic_cast<const Type::SelectionBase*>(type);
        if (selection_type != NULL ) return make_storage(p, selection_type );

        // default error
    }

}



StorageBase * JSONToStorage::make_storage(JSONPath &p, const Type::Record *record)
{
    if (p.head()->type() == json_spirit::obj_type) {
        // check individual keys
        for( Type::Record::KeyIter it= record->begin(); it != record->end(); ++it) {

        }
    }
    // possibly construction of reduced record

}



StorageBase * JSONToStorage::make_storage(JSONPath &p, const Type::AbstractRecord *abstr_rec)
{}



StorageBase * JSONToStorage::make_storage(JSONPath &p, const Type::Array *array)
{
    if (p.head()->type() == json_spirit::array_type) {
        const json_spirit::mArray & j_array = p.head()->get_array();
        if ( ! array->match_size( j_array.size() ) ) {
            THROW( ExcInputError() << EI_Specification("Do not fit into size limits of the Array.") << EI_ErrorAddress(p) << EI_InputType(array) );
        }
    } else {
        THROW( ExcInputError() << EI_Specification("Wrong type, has to be Array.") << EI_ErrorAddress(p) << EI_InputType(array) );
    }
}



StorageBase * JSONToStorage::make_storage(JSONPath &p, const Type::SelectionBase *selection)
{
    if (p.head()->type() == json_spirit::str_type) {
        try {
            int value = selection->name_to_int( p.head()->get_str() );
            return new StorageInt( value );
        } catch (Type::SelectionBase::ExcSelectionKeyNotFound &exc) {
            THROW( ExcInputError() << EI_Specification("Wrong value of the Selection.") << EI_ErrorAddress(p) << EI_InputType(selection) );
        }
    }
    THROW( ExcInputError() << EI_Specification("Wrong type, value should be String (key of Selection).\n") << EI_ErrorAddress(p) << EI_InputType(selection) );
    return NULL;
}



StorageBase * JSONToStorage::make_storage(JSONPath &p, const Type::Bool *bool_type)
{
    if (p.head()->type() == json_spirit::bool_type) {
        return new StorageBool( p.head()->get_bool() );
    } else {
        xprintf(UsrErr,"The input value at address '%s' has to be of type Bool.\n", p.str().c_str());
    }
    return NULL;
}



StorageBase * JSONToStorage::make_storage(JSONPath &p, const Type::Integer *int_type)
{
    if (p.head()->type() == json_spirit::int_type) {
        int value = p.head()->get_int();
        if (int_type->match(value)) {
            return new StorageInt( value );
        } else {
            xprintf(UsrErr,"The input value at address '%s' has to be of type Int.\n", p.str().c_str());
        }

    } else {
        xprintf(UsrErr,"The input value at address '%s' out of bounds.\n", p.str().c_str());
    }
    return NULL;
}



StorageBase * JSONToStorage::make_storage(JSONPath &p, const Type::Double *double_type)
{
    if (p.head()->type() == json_spirit::real_type) {
        double value = p.head()->get_real();
        double_type->match(value);
        return new StorageDouble( value );
    } else {
        xprintf(UsrErr,"The input value at address '%s' has to be of type Double.\n", p.str().c_str());
    }
    return NULL;
}



StorageBase * JSONToStorage::make_storage(JSONPath &p, const Type::String *string_type)
{
    if (p.head()->type() == json_spirit::str_type) {
        string value = p.head()->get_str();
        //double_type->match(value);        // possible parsing and modifications of special strings
        return new StorageString( value );
    } else {
        xprintf(UsrErr,"The input value at address '%s' has to be of type String.\n", p.str().c_str());
    }
    return NULL;
}



} // namespace Input
