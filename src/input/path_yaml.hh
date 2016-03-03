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
 * @file    path_yaml.hh
 * @brief   
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
	/// Definition of YAML-cpp node.
    typedef YAML::Node Node;

    /// Constructor.
    PathYAML(std::istream &in);

    /**
     * @brief Destructor.
     *
     * Clean nodes_.
     */
    ~PathYAML() override;

    /**
     * @brief Returns level of actual path.
     *
     * Root has level == 0.
     */
    inline int level() const
    { return nodes_.size() - 1; }

    /**
     * @brief Dive into yaml-cpp hierarchy.
     *
     * Implements @p PathBase::down(unsigned int)
     */
    bool down(unsigned int index) override;

    /**
     * @brief Dive into yaml-cpp hierarchy.
     *
     * Implements @p PathBase::down(const std::string&)
     */
    bool down(const std::string& key) override;

    /// Return one level up in the hierarchy.
    void up() override;

    // These methods are derived from PathBase
    bool is_null_type() const override;                               ///< Implements @p PathBase::is_null_type
    bool get_bool_value() const override;                             ///< Implements @p PathBase::get_bool_value
    std::int64_t get_int_value() const override;                      ///< Implements @p PathBase::get_int_value
    double get_double_value() const override;                         ///< Implements @p PathBase::get_double_value
    std::string get_string_value() const override;                    ///< Implements @p PathBase::get_string_value
    unsigned int get_node_type_index() const override;                ///< Implements @p PathBase::get_node_type_index
    bool get_record_key_set(std::set<std::string> &) const override;  ///< Implements @p PathBase::get_record_key_set
    int get_array_size() const override;                              ///< Implements @p PathBase::get_array_size
    bool is_record_type() const override;                             ///< Implements @p PathBase::is_record_type
    bool is_array_type() const override;                              ///< Implements @p PathBase::is_array_type
    PathYAML * clone() const override;                                ///< Implements @p PathBase::clone

    /// Implements reading of reference keys, and check of cyclic references.
    PathBase * find_ref_node() override;

    /// Implements @p PathBase::get_descendant_name
    std::string get_descendant_name() const override;

protected:


    /// Pointer to YAML Value object at current path.
    inline const Node &head() const
    { return nodes_.back(); }

    // Pointers to all nodes from the root up to the current path.
    std::vector<Node> nodes_;
};


/**
 * @brief Output operator for PathYAML.
 *
 * Mainly for debugging purposes and error messages.
 */
std::ostream& operator<<(std::ostream& stream, const PathYAML& path);



} // namespace Input



#endif /* PATH_YAML_HH_ */
