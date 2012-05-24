/*
 * json_to_storage.hh
 *
 *  Created on: May 7, 2012
 *      Author: jb
 *
 * TODO:
 *  EI_InputType - passing pointer is not save. The object can be deleted during stack unfolding.
 *  Error_info can not be used for abstract types (of course), the only save way is to
 *  use  boost::shared_ptr, so make the copy our self by make_shared and then only pass the pointer.
 *
 *  These problems can be eliminated, when all Type instances are static.
 *
 * TODO:
 * - check cyclic references, drop const for json_spirit pointers and modify REF keys
 *   when dereferenced and modify it back when we return.
 *   (e.g.  add a new entry into map
 * TODO:
 *   Find deleting of a null pointer in json_spirit library. Possibly implement
 *   output of true stack for GNU.
 */

#ifndef JSON_TO_STORAGE_HH_
#define JSON_TO_STORAGE_HH_



#include "input/input_type.hh"
#include "input/input_interface.hh"
#include "input/storage.hh"
#include "json_spirit/json_spirit.h"


namespace Input {



/**
 * This class represents an address in the JSON file.
 */
class JSONPath {
public:
    /**
     * Thrown if a reference in the input file
     */

    TYPEDEF_ERR_INFO(EI_ErrorAddress, JSONPath);
    TYPEDEF_ERR_INFO(EI_RefAddress, JSONPath);
    TYPEDEF_ERR_INFO(EI_RefStr, const string);
    TYPEDEF_ERR_INFO(EI_Specification, const string);
    DECLARE_EXCEPTION(ExcRefOfWrongType,
            << "Reference at address "
            << EI_ErrorAddress::qval << " has wrong type, should by string.");
    DECLARE_EXCEPTION(ExcReferenceNotFound,
            << "Reference {REF=\"" << EI_RefStr::val << "\"} at address " << EI_RefAddress::qval << " not found.\n"
            << "failed to follow at address: " << EI_ErrorAddress::qval << " because " << EI_Specification::val);



    typedef json_spirit::mValue Node;

    JSONPath(const Node& root_node);

    /**
     * Dive into json_spirit hierarchy. Store current path and returns pointer to new json_spirit node.
     * If the json_spirit type do not match returns NULL.
     */
    const Node * down(unsigned int index);
    const Node * down(const string& key);

    /**
     * Return one level up in the hierrarchy.
     */
    void up();

    /**
     * Move to root node.
     */
    void go_to_root();

    /**
     * Pointer to JSON Value object at current path.
     */
    inline const Node * head() const
    { return nodes_.back(); }

    /**
     * Returns level of actual path. Root has level == 0.
     */
    inline int level() const
    { return nodes_.size() - 1; }

    /**
     * Check if current head node is a JSON Object containing one key REF of type string.
     * If yes, returns the string through reference @p ref_address.
     */
    bool get_ref_from_head(string & ref_address);

    /**
     * Looks for key 'TYPE' and possibly return its value Node.
     */
    const Node * get_abstract_type_from_head();

    /**
     * Creates a new JSONPath object given by  address string possibly relative to the current
     * path.
     */
    JSONPath find_ref_node(const string& ref_address);
    /**
     * Output to the given stream.
     */
    void output(ostream &stream) const;
    /**
     * Returns string address of current position.
     */
    string str();

private:
    /**
     * One level of the @p path_ is either index (nonnegative int) in array or string key in a json object.
     * For the first type we save index into first part of the pair and empty string to the second.
     * For the later type of level, we save -1 for index and the key into the secodn part of the pair.
     */
    vector< pair<int, string> > path_;
    vector<const Node *> nodes_;

};

std::ostream& operator<<(std::ostream& stream, const JSONPath& path);


/**
 *
 */
class JSONToStorage {
public:
    /**
     * Exceptions.
     */
    TYPEDEF_ERR_INFO(EI_InputType, Type::TypeBase const *);
    TYPEDEF_ERR_INFO(EI_File, const string);
    TYPEDEF_ERR_INFO(EI_Specification, const string);
    TYPEDEF_ERR_INFO( EI_ErrorAddress, JSONPath);
    DECLARE_EXCEPTION( ExcInputError, << "Error in input file: " << EI_File::val << " at address: " << EI_ErrorAddress::val <<"\n"
                                            << EI_Specification::val << "\n"
                                            << "Expected type:\n" << *EI_InputType::ref(_exc) );

    TYPEDEF_ERR_INFO( EI_JSONLine, unsigned int);
    TYPEDEF_ERR_INFO( EI_JSONColumn, unsigned int);
    TYPEDEF_ERR_INFO( EI_JSONReason, string);
    DECLARE_EXCEPTION( ExcNotJSONFormat, << "Not valid JSON file " << EI_File::qval << ". Error at line "
            << EI_JSONLine::val << " : col " << EI_JSONColumn::val
            << " ; reason: " << EI_JSONReason::val << "\n" );

    JSONToStorage();

    void read_stream(istream &in, const Type::TypeBase &root_type);

    template <class T>
    T get_root_interface() const;

    const StorageBase *get_storage()
    { return storage_;}

protected:
    /**
     * Check correctness of the input given by json_spirit node at head() of JSONPath @p p
     * against type specification @p type. Die on input error (and return NULL).
     * For correct input store values into storage and return pointer to root of created storage tree.
     */
    StorageBase * make_storage(JSONPath &p, const Type::TypeBase *type);

    StorageBase * make_storage(JSONPath &p, const Type::Record *record);
    StorageBase * make_storage(JSONPath &p, const Type::AbstractRecord *abstr_rec);
    StorageBase * make_storage(JSONPath &p, const Type::Array *array);
    StorageBase * make_storage(JSONPath &p, const Type::SelectionBase *selection);
    StorageBase * make_storage(JSONPath &p, const Type::Bool *bool_type);
    StorageBase * make_storage(JSONPath &p, const Type::Integer *int_type);
    StorageBase * make_storage(JSONPath &p, const Type::Double *double_type);
    StorageBase * make_storage(JSONPath &p, const Type::String *string_type);



    StorageBase *storage_;
    StorageArray *envelope; // helper envelope for get_root_interface
    const Type::TypeBase *root_type_;
};







/********************************************88
 * Implementation
 */

template <class T>
T JSONToStorage::get_root_interface() const
{
    return *(Iterator<T>( *root_type_, envelope, 0));
}




} // namespace Input

class Some {

};

#endif /* JSON_TO_STORAGE_HH_ */
