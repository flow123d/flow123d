/*
 * path_base.hh
 *
 *  Created on: May 7, 2012
 *      Author: jb
 */

#ifndef PATH_BASE_HH_
#define PATH_BASE_HH_


#include "input/accessors.hh"


namespace Input {
using namespace std;


/**
 * @brief Base abstract class used by ReaderToStorage class to iterate over the input tree.
 *
 * Currently this class has two descendants
 *  - PathJSON: work with JSON input tree
 *  - PathYAML: work with YAML input tree
 */
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

    /// Get short string description of node type, method is used for printout of messages
    std::string get_node_type(unsigned int type_idx) const;

    /// Get index of head type, value corresponds with order in @p json_type_names vector
    virtual unsigned int get_node_type_index() const =0;

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

protected:
    PathBase();

    /**
     * One level of the @p path_ is either index (nonnegative int) in array or string key in a json object.
     * For the first type we save index into first part of the pair and empty string to the second.
     * For the later type of level, we save -1 for index and the key into the secodn part of the pair.
     */
    vector< pair<int, string> > path_;

    /**
     * Names of all possible node types in parsed JSON tree provided by JSON Spirit library.
     * Initialized in constructor.
     *
     */
    vector<string> json_type_names;

};



} // namespace Input



#endif /* PATH_BASE_HH_ */
