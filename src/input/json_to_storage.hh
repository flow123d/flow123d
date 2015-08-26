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
#include "yaml-cpp/yaml.h"
#include "input/input_type.hh"

#include "input/accessors.hh"
#include "input/storage.hh"


namespace Input {



/// Possible formats of input files.
typedef enum  {
    format_JSON,
    format_YAML
} FileFormat;



class PathBase {
public:

    /**
     * Thrown if a reference in the input file
     */
    TYPEDEF_ERR_INFO(EI_ErrorAddress, string);
    TYPEDEF_ERR_INFO(EI_RefAddress, string);
    TYPEDEF_ERR_INFO(EI_JsonFile, const string);
    TYPEDEF_ERR_INFO(EI_RefStr, const string);
    TYPEDEF_ERR_INFO(EI_Specification, const string);
    DECLARE_INPUT_EXCEPTION(ExcRefOfWrongType,
            << "Reference at address "
            << EI_ErrorAddress::qval << " has wrong type, should by string.");
    DECLARE_INPUT_EXCEPTION(ExcReferenceNotFound,
            << "Error in input file: " << EI_JsonFile::qval << "\nReference {REF=\"" << EI_RefStr::val << "\"} at address " << EI_RefAddress::qval << " not found.\n"
            << "failed to follow at address: " << EI_ErrorAddress::qval << " because " << EI_Specification::val);



    /**
     * Require virtual destructor also for child classes.
     */
    virtual ~PathBase() {};

    /**
     * Returns level of actual path. Root has level == 0.
     */
	virtual int level() const =0;

    /**
     * Check if current head node is containing one key REF of type string.
     *
     * If yes, creates a new path object given by address string possibly relative to the current
     * path. In other else return NULL.
     *
     * This method has the meaning only for JSON. For YAML (YAML has native references) return
     * always NULL.
     */
	virtual PathBase * find_ref_node() =0;

	/**
	 * Create copy of derived class.
	 */
	virtual PathBase * clone() const =0;

    /**
     * Output to the given stream.
     */
    void output(ostream &stream) const;

    /// Check if type of head node is null
    virtual bool is_null_type() const =0;

    /// Get boolean value of head node or throw exception
    virtual bool get_bool_value() const =0;

    /// Get integer value of head node or throw exception
    virtual std::int64_t get_int_value() const =0;

    /// Get double value of head node or throw exception
    virtual double get_double_value() const =0;

    /// Get string value of head node or throw exception
    virtual std::string get_string_value() const =0;

    /// Get short string description of head type, method is used for printout of messages
    virtual std::string get_node_type() const =0;

    /// Get set of keys of head type record, if head type is not record return false
    virtual bool get_record_key_set(std::set<std::string> &) const =0;

    /// Get size of array (sequence type), if object is not array return -1
    virtual int get_array_size() const =0;

    /// Check if type of head node is record
    virtual bool is_record_type() const =0;

    /// Check if type of head node is array
    virtual bool is_array_type() const =0;

    /**
     * Dive into json_spirit hierarchy. Store current path and returns true if pointer to new json_spirit node is not NULL.
     */
    virtual bool down(unsigned int index) =0;
    virtual bool down(const string& key) =0;

    /**
     * Return one level up in the hierarchy.
     */
    virtual void up() =0;

    /**
     * Move to root node.
     */
    void go_to_root();

    /**
     * Returns string address of current position.
     */
    string as_string() const;

    /**
     * Gets name of descendant of Abstract:
     * - for JSON returns value of TYPE key
     * - for YAML returns value of tag
     *
     * If descendant name is not found returns empty string.
     */
    virtual std::string get_descendant_name() const =0;

    /// Returns name of sequence type (e. g. "JSON array")
    virtual std::string sequence_name() const =0;

    /// Returns name of map type (e. g. "JSON object")
    virtual std::string map_name() const =0;

protected:
    PathBase();

    /**
     * One level of the @p path_ is either index (nonnegative int) in array or string key in a json object.
     * For the first type we save index into first part of the pair and empty string to the second.
     * For the later type of level, we save -1 for index and the key into the secodn part of the pair.
     */
    vector< pair<int, string> > path_;

};



/**
 * @brief Class used by JSONToStorage class to iterate over the JSON tree provided by json_spirit library.
 *
 * This class keeps whole path from the root of the JSON tree to the current node. We store nodes along path in \p nodes_
 * and address of the node in \p path_.
 *
 * The class also contains methods for processing of special keys 'REF' and 'TYPE'. The reference is record with only one key
 * 'REF' with a string value that contains address of the reference. The string with the address is extracted and provided by
 * method \p JSONtoStorage::find_ref_node.
 */
class PathJSON : public PathBase {
public:
    typedef json_spirit::mValue Node;

    PathJSON(istream &in);

    /**
     * Dive into json_spirit hierarchy. Store current path and returns true if pointer to new json_spirit node is not NULL.
     */
    bool down(unsigned int index) override;
    bool down(const string& key) override;

    /**
     * Return one level up in the hierarchy.
     */
    void up() override;

    /**
     * Returns level of actual path. Root has level == 0.
     */
    inline int level() const
    { return nodes_.size() - 1; }

    // These methods are derived from PathBase
    bool is_null_type() const override;
    bool get_bool_value() const override;
    std::int64_t get_int_value() const override;
    double get_double_value() const override;
    std::string get_string_value() const override;
    std::string get_node_type() const override;
    bool get_record_key_set(std::set<std::string> &) const override;
    int get_array_size() const override;
    bool is_record_type() const override;
    bool is_array_type() const override;
    PathJSON * clone() const override;

    inline std::string sequence_name() const override {
    	return "JSON array";
    }

    inline std::string map_name() const override {
    	return "JSON object";
    }

    PathBase * find_ref_node();

    /**
     * Put address of actual reference to previous_references_ set
     */
    void remember_reference();

    std::string get_descendant_name() const;

protected:

    /**
     * Default constructor.
     * Provides common initialization for public constructors.
     */
    PathJSON();

    /**
     * Pointer to JSON Value object at current path.
     */
    inline const Node * head() const
    { return nodes_.back(); }

    std::set<string> previous_references_;

    vector<const Node *> nodes_;

    /**
     * Names of all possible node types in parsed JSON tree provided by JSON Spirit library.
     * Initialized in constructor.
     *
     */
    vector<string> json_type_names;

};

/**
 * Output operator for PathJSON. Mainly for debugging purposes and error messages.
 */
std::ostream& operator<<(std::ostream& stream, const PathJSON& path);


class PathYAML : public PathBase {
public:
    typedef YAML::Node Node;

    PathYAML(istream &in);

    /**
     * Returns level of actual path. Root has level == 0.
     */
    inline int level() const
    { return nodes_.size() - 1; }

    /**
     * Dive into yaml-cpp hierarchy. Store current path and returns true if pointer to new yaml node is not NULL.
     */
    bool down(unsigned int index) override;
    bool down(const string& key) override;

    /**
     * Return one level up in the hierarchy.
     */
    void up() override;

    // These methods are derived from PathBase
    bool is_null_type() const override;
    bool get_bool_value() const override;
    std::int64_t get_int_value() const override;
    double get_double_value() const override;
    std::string get_string_value() const override;
    std::string get_node_type() const override;
    bool get_record_key_set(std::set<std::string> &) const override;
    int get_array_size() const override;
    bool is_record_type() const override;
    bool is_array_type() const override;
    PathYAML * clone() const override;

    inline std::string sequence_name() const override {
    	return "YAML sequence";
    }

    inline std::string map_name() const override {
    	return "YAML map";
    }

    PathBase * find_ref_node();

    std::string get_descendant_name() const;

protected:
    /**
     * Pointer to YAML Value object at current path.
     */
    inline const Node * head() const
    { return nodes_.back(); }

    vector<const Node *> nodes_;
};


/**
 * Output operator for PathYAML. Mainly for debugging purposes and error messages.
 */
std::ostream& operator<<(std::ostream& stream, const PathYAML& path);


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
    /// General exception during conversion from JSON tree to storage.
    TYPEDEF_ERR_INFO(EI_InputType, string );
    TYPEDEF_ERR_INFO(EI_File, const string);
    TYPEDEF_ERR_INFO(EI_Specification, const string);
    TYPEDEF_ERR_INFO(EI_Format, const string);
    TYPEDEF_ERR_INFO(EI_JSON_Type, const string);
    TYPEDEF_ERR_INFO( EI_ErrorAddress, string);
    DECLARE_INPUT_EXCEPTION( ExcInputError, << "Error in input file: " << EI_File::qval << " at address: " << EI_ErrorAddress::qval << "\n"
                                            << EI_Specification::val << "\n"
                                            << EI_Format::val << " type: " << EI_JSON_Type::qval << "\n"
                                            << "Expected type:\n" << EI_InputType::val );


    TYPEDEF_ERR_INFO( EI_JSONLine, unsigned int);
    TYPEDEF_ERR_INFO( EI_JSONColumn, unsigned int);
    TYPEDEF_ERR_INFO( EI_JSONReason, string);
    DECLARE_INPUT_EXCEPTION( ExcNotJSONFormat, << "Not valid JSON file " << EI_File::qval << ". Error at line "
            << EI_JSONLine::val << " : col " << EI_JSONColumn::val
            << " ; reason: " << EI_JSONReason::val << "\n" );

    TYPEDEF_ERR_INFO( EI_InputErrorMessage, const string);
    TYPEDEF_ERR_INFO( EI_RecordName, const string);
    DECLARE_INPUT_EXCEPTION( ExcAutomaticConversionError, << "Error during automatic conversion of "
    		<< EI_RecordName::val << " record.\n " << EI_InputErrorMessage::val << "\n" );


    /**
     * Read a storage from input stream. Parameter @p root_type
     * provides input type tree declaration. See @p read_from_stream for details.
     */
    JSONToStorage(const FilePath &in_file, const Type::TypeBase &root_type);

    /**
     * Read a storage from string (e.g. complex default value).
     */
    JSONToStorage( const string &default_str, const Type::TypeBase &root_type, FileFormat format);

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
    void read_stream(istream &in, const Type::TypeBase &root_type, FileFormat format);

    /**
     * Getter for root of the storage tree.
     */
    const StorageBase *get_storage()
    { return storage_;}


    /**
     * Check correctness of the input given by json_spirit node at head() of PathJSON @p p
     * against type specification @p type. Die on input error (and return NULL).
     * For correct input, creates the storage tree and returns pointer to its root node.
     */
    StorageBase * make_storage(PathBase &p, const Type::TypeBase *type);

    StorageBase * make_storage(PathBase &p, const Type::Record *record);
    StorageBase * make_storage(PathBase &p, const Type::AbstractRecord *abstr_rec);
    StorageBase * make_storage(PathBase &p, const Type::Array *array);

    StorageBase * make_selection_storage_without_catch(PathBase &p, const Type::Selection *selection);
    StorageBase * make_storage(PathBase &p, const Type::Selection *selection);
    StorageBase * make_storage(PathBase &p, const Type::Bool *bool_type);
    StorageBase * make_storage(PathBase &p, const Type::Integer *int_type);
    StorageBase * make_storage(PathBase &p, const Type::Double *double_type);
    StorageBase * make_storage(PathBase &p, const Type::String *string_type);

    StorageBase * record_automatic_conversion(PathBase &p, const Type::Record *record);
    StorageBase * abstract_rec_automatic_conversion(PathBase &p, const Type::AbstractRecord *abstr_rec);

    /**
     * Dispatch according to @p type and create corresponding storage from the given string.
     */
    StorageBase * make_storage_from_default( const string &dflt_str, const Type::TypeBase *type);


    /// Storage of the read and checked input data
    StorageBase *storage_;

    /// Root of the declaration tree of the data in the storage.
    const Type::TypeBase *root_type_;

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
