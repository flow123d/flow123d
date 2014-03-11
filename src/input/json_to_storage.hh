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


#include <sstream>

#include "json_spirit/json_spirit.h"
#include "input/input_type.hh"

#include "input/accessors.hh"
#include "input/storage.hh"


namespace Input {



/**
 * @brief Class used by JSONToStorage class to iterate over the JSON tree provided by json_spirit library.
 *
 * This class keeps whole path from the root of the JSON tree to the current node. We store nodes along path in \p nodes_
 * and address of the node in \p path_.
 *
 * The class also contains methods for processing of special keys 'REF' and 'TYPE'. The reference is record with only one key
 * 'REF' with a string value that contains address of the reference. The string with the address is extracted by \p JSONToStorage::get_ref_from_head
 * then the JSONPath corresponding to the address is provided by method \p JSONtoStorage::find_ref_node.
 */
class JSONPath {
public:
    /**
     * Thrown if a reference in the input file
     */

    TYPEDEF_ERR_INFO(EI_ErrorAddress, JSONPath);
    TYPEDEF_ERR_INFO(EI_RefAddress, JSONPath);
    TYPEDEF_ERR_INFO(EI_JsonFile, const string);
    TYPEDEF_ERR_INFO(EI_RefStr, const string);
    TYPEDEF_ERR_INFO(EI_Specification, const string);
    DECLARE_INPUT_EXCEPTION(ExcRefOfWrongType,
            << "Reference at address "
            << EI_ErrorAddress::qval << " has wrong type, should by string.");
    DECLARE_INPUT_EXCEPTION(ExcReferenceNotFound,
            << "Error in input file: " << EI_JsonFile::qval << "\nReference {REF=\"" << EI_RefStr::val << "\"} at address " << EI_RefAddress::qval << " not found.\n"
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

    /**
     * Put actual address to previous_references_ set
     */
    void put_address();

private:
    /**
     * One level of the @p path_ is either index (nonnegative int) in array or string key in a json object.
     * For the first type we save index into first part of the pair and empty string to the second.
     * For the later type of level, we save -1 for index and the key into the secodn part of the pair.
     */
    vector< pair<int, string> > path_;
    vector<const Node *> nodes_;
    std::set<string> previous_references_;


};

/**
 * Output operator for JSONPath. Mainly for debugging purposes and error messages.
 */
std::ostream& operator<<(std::ostream& stream, const JSONPath& path);


/**
 *  @brief Reader for (slightly) modified JSON files.
 *
 *  This class implements a reader of modified JSON file format. The modifications include
 *  shell-like comments (using hash '#' character), this is implemented in comment_filter.hh,  optional quoting of
 *  keys in JSON objects that do not contain spaces, and possibility to use '=' instead of ':'. So you can write:
 *  @code
 *    { key1="text", key2=2, "key 3"=4 }
 *  @endcode
 *  Note, however, that our input interface allows only C identifiers for keys. The reader use json_spirit library
 *  (based on Spirit parser from Boost) with slightly modified grammar.
 *
 *  The input file is at first read and parsed by json_spirit. Then JSONToStorage pass through tree with parsed data along
 *  with passing through declaration tree. The input data are check against declaration and stored in the Storage tree.
 *
 *  Accessor to the root record is provided by JSONToStorage::get_root_interface<T> method template.
 *
 *  @ingroup input
 */
class JSONToStorage {
public:
    /*
     * Exceptions.
     */
    // unfortunately following is not safe:
    // TYPEDEF_ERR_INFO(EI_InputType, Type::TypeBase const *);

    /// General exception during conversion from JSON tree to storage.
    TYPEDEF_ERR_INFO(EI_InputType, string );
    TYPEDEF_ERR_INFO(EI_File, const string);
    TYPEDEF_ERR_INFO(EI_Specification, const string);
    TYPEDEF_ERR_INFO(EI_JSON_Type, const string);
    TYPEDEF_ERR_INFO( EI_ErrorAddress, JSONPath);
    DECLARE_INPUT_EXCEPTION( ExcInputError, << "Error in input file: " << EI_File::qval << " at address: " << EI_ErrorAddress::qval <<"\n"
                                            << EI_Specification::val << "\n"
                                            << "JSON type: " << EI_JSON_Type::qval << "\n"
                                            << "Expected type:\n" << EI_InputType::val );


    TYPEDEF_ERR_INFO( EI_JSONLine, unsigned int);
    TYPEDEF_ERR_INFO( EI_JSONColumn, unsigned int);
    TYPEDEF_ERR_INFO( EI_JSONReason, string);
    DECLARE_INPUT_EXCEPTION( ExcNotJSONFormat, << "Not valid JSON file " << EI_File::qval << ". Error at line "
            << EI_JSONLine::val << " : col " << EI_JSONColumn::val
            << " ; reason: " << EI_JSONReason::val << "\n" );


    /**
     * Read a storage from input stream. Parameter @p root_type
     * provides input type tree declaration. See @p read_from_stream for details.
     */
    JSONToStorage(istream &in, const Type::TypeBase &root_type);

    /**
     * Read a storage from string (e.g. complex default value).
     */
    JSONToStorage( const string &default_str, const Type::TypeBase &root_type);

    /**
     * Returns the root accessor. The template type \p T should correspond
     * to the kind of the input type at root of the declaration tree.
     */
    template <class T>
    T get_root_interface() const;


protected:

    /**
     * Default constructor.
     * Provides common initialization for public constructors.
     */
    JSONToStorage();

    /**
     * This method actually reads the given stream \p in, checks the data just read against the declaration tree given by \p root_type, and
     * store the data into private storage tree using \p StorageBase classes.
     */
    void read_stream(istream &in, const Type::TypeBase &root_type);

    /**
     * Getter for root of the storage tree.
     */
    const StorageBase *get_storage()
    { return storage_;}


    /**
     * Check correctness of the input given by json_spirit node at head() of JSONPath @p p
     * against type specification @p type. Die on input error (and return NULL).
     * For correct input, creates the storage tree and returns pointer to its root node.
     */
    StorageBase * make_storage(JSONPath &p, const Type::TypeBase *type);

    StorageBase * make_storage(JSONPath &p, const Type::Record *record);
    StorageBase * make_storage(JSONPath &p, const Type::AbstractRecord *abstr_rec);
    StorageBase * make_storage(JSONPath &p, const Type::Array *array);

    StorageBase * make_selection_storage_without_catch(JSONPath &p, const Type::Selection *selection);
    StorageBase * make_storage(JSONPath &p, const Type::Selection *selection);
    StorageBase * make_storage(JSONPath &p, const Type::Bool *bool_type);
    StorageBase * make_storage(JSONPath &p, const Type::Integer *int_type);
    StorageBase * make_storage(JSONPath &p, const Type::Double *double_type);
    StorageBase * make_storage(JSONPath &p, const Type::String *string_type);

    /**
     * Dispatch according to @p type and create corresponding storage from the given string.
     */
    StorageBase * make_storage_from_default( const string &dflt_str, const Type::TypeBase *type);



    /// helper envelope for get_root_interface
    //StorageArray *envelope;

    /// Storage of the read and checked input data
    StorageBase *storage_;

    /// Root of the declaration tree of the data in the storage.
    const Type::TypeBase *root_type_;

    /**
     * Names of all possible node types in parsed JSON tree provided by JSON Spirit library.
     * Initialized in constructor.
     *
     */
    vector<string> json_type_names;

    /**
     * List of paths specifications of the keys that wasn't read by the JSON reader.
     */
    //vector<string> ignored_keys;

};







/********************************************88
 * Implementation
 */

template <class T>
T JSONToStorage::get_root_interface() const
{
	ASSERT(storage_, "NULL pointer to storage !!! \n");

	Address addr(storage_, root_type_);
	// try to create an iterator just to check type
	Iterator<T>( *root_type_, addr, 0);

	auto tmp_root_type = static_cast<const typename T::InputType &>(*root_type_);
    return T( addr, tmp_root_type );
}




} // namespace Input



#endif /* JSON_TO_STORAGE_HH_ */
