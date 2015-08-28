/*
 * path_json.hh
 *
 *  Created on: May 7, 2012
 *      Author: jb
 */

#ifndef PATH_JSON_HH_
#define PATH_JSON_HH_


#include "input/path_base.hh"
#include "json_spirit/json_spirit.h"


namespace Input {



/**
 * @brief Class used by ReaderToStorage class to iterate over the JSON tree provided by json_spirit library.
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
    unsigned int get_node_type_index() const override;
    bool get_record_key_set(std::set<std::string> &) const override;
    int get_array_size() const override;
    bool is_record_type() const override;
    bool is_array_type() const override;
    PathJSON * clone() const override;

    PathBase * find_ref_node() override;

    /**
     * Put address of actual reference to previous_references_ set
     */
    void remember_reference();

    std::string get_descendant_name() const override;

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

};

/**
 * Output operator for PathJSON. Mainly for debugging purposes and error messages.
 */
std::ostream& operator<<(std::ostream& stream, const PathJSON& path);



} // namespace Input



#endif /* PATH_JSON_HH_ */
