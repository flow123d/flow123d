/*
 * path_yaml.hh
 *
 *  Created on: May 7, 2012
 *      Author: jb
 */

#ifndef PATH_YAML_HH_
#define PATH_YAML_HH_


#include "yaml-cpp/yaml.h"
#include "input/path_base.hh"


namespace Input {



/**
 * @brief Class used by ReaderToStorage class to iterate over the YAML tree provided by yaml-cpp library.
 *
 * This class keeps whole path from the root of the YAML tree to the current node. We store nodes along path in \p nodes_
 * and address of the node in \p path_.
 *
 * The class also contains methods for processing of special tags for 'TYPE' key. Class doesn't need special methods
 * for work with references. YAML tree used own native references.
 */
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
    unsigned int get_node_type_index() const override;
    bool get_record_key_set(std::set<std::string> &) const override;
    int get_array_size() const override;
    bool is_record_type() const override;
    bool is_array_type() const override;
    PathYAML * clone() const override;

    PathBase * find_ref_node() override;

    std::string get_descendant_name() const override;

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



} // namespace Input



#endif /* PATH_YAML_HH_ */
