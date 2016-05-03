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
 * @file    path_base.hh
 * @brief   
 */

#ifndef PATH_BASE_HH_
#define PATH_BASE_HH_

#include <utility>
#include <set>
#include <string>
#include <vector>

#include "input/input_exception.hh"


namespace Input {


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
    TYPEDEF_ERR_INFO(EI_ErrorAddress, std::string);
    TYPEDEF_ERR_INFO(EI_RefAddress, std::string);
    TYPEDEF_ERR_INFO(EI_JsonFile, std::string);
    TYPEDEF_ERR_INFO(EI_RefStr, std::string);
    TYPEDEF_ERR_INFO(EI_Specification, std::string);
    DECLARE_INPUT_EXCEPTION(ExcRefOfWrongType,
            << "Reference at address "
            << EI_ErrorAddress::qval << " has wrong type, should by string.");
    DECLARE_INPUT_EXCEPTION(ExcReferenceNotFound,
            << "Error in input file: " << EI_JsonFile::qval << "\nReference {REF=\"" << EI_RefStr::val << "\"} at address " << EI_RefAddress::qval << " not found.\n"
            << "failed to follow at address: " << EI_ErrorAddress::qval << " because " << EI_Specification::val);

    /// Must have virtual destructor to call the right one form child.
    virtual ~PathBase() {};


    /**
     * @brief Returns level of actual path.
     *
     * Root has level == 0.
     */
	virtual int level() const =0;

    /**
     * @brief Check if current head node is containing one key REF of type string.
     *
     * If yes, creates a new path object given by address string possibly relative to the current
     * path. In other else return NULL.
     *
     * This method has the meaning only for JSON. For YAML (YAML has native references) return
     * always NULL.
     */
	virtual PathBase * find_ref_node() =0;

	/// Create copy of derived class.
	virtual PathBase * clone() const =0;

    /// Output to the given stream.
    void output(std::ostream &stream) const;

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
     * @brief Dive one level down into path hierarchy.
     *
     * Store current path and returns true if pointer to new node is not NULL.
     */
    virtual bool down(unsigned int index) =0;

    /**
     * @brief Dive one level down into path hierarchy.
     *
     * Store current path and returns true if pointer to new node is not NULL.
     */
    virtual bool down(const std::string& key) =0;

    /// Return one level up in the hierarchy.
    virtual void up() =0;

    /// Move to root node.
    void go_to_root();

    /// Returns string address of current position.
    std::string as_string() const;

    /**
     * @brief Gets name of descendant of Abstract.
     *
     * - for JSON returns value of TYPE key
     * - for YAML returns value of tag
     *
     * If descendant name is not found returns empty string.
     */
    virtual std::string get_descendant_name() const =0;

protected:
    /// Forbid default constructor.
    PathBase();

    /**
     * @brief One level of the @p path_ is either index (nonnegative int) in array or string key in a json object.
     *
     * For the first type we save index into first part of the pair and empty string to the second.
     * For the later type of level, we save -1 for index and the key into the secodn part of the pair.
     */
    std::vector< std::pair<int, std::string> > path_;

    /**
     * @brief Names of all possible node types in parsed input tree.
     *
     * Names are provided by JSON Spirit or YAML-cpp library.
     * Initialized in constructor.
     *
     */
    std::vector<std::string> json_type_names;

};



} // namespace Input



#endif /* PATH_BASE_HH_ */
