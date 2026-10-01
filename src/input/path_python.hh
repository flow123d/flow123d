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
 * @file    path_python.hh
 * @brief
 */

#ifndef PATH_PYTHON_HH_
#define PATH_PYTHON_HH_

#pragma once

#include <memory>
#include <stdint.h>                               // for int64_t
#include <iosfwd>                                 // for ostream
#include <set>                                    // for set
#include <string>                                 // for string
#include <vector>                                 // for vector
#include <pybind11/pybind11.h>
#include "input/path_base.hh"

namespace py = pybind11;

namespace Input {

// Pybind11 needs set visibility to hidden (see https://pybind11.readthedocs.io/en/stable/faq.html).
#pragma GCC visibility push(hidden)

/**
 * @brief Class used by ReaderToStorage class to iterate over the tree provided by Python dict and list objects.
 *
 * This class keeps whole path from the root of the Python objects tree to the current node. We store nodes along path
 * in \p nodes_ and address of the node in \p path_.
 *
 * The class also contains methods for processing of special keys 'REF' and 'TYPE'. The reference is record with only one key
 * 'REF' with a string value that contains address of the reference. The string with the address is extracted and provided by
 * method \p JSONtoStorage::find_ref_node.
 */
class PathPython : public PathBase {
public:

    enum class PythonNodeType {
        dict = 0,
        list,
        string,
        boolean,
        integer,
        real,
        none,
        other,
        undefined
    };

    explicit PathPython(const py::dict &root);

    ~PathPython() override;

    bool down(unsigned int index) override;
    bool down(const std::string &key, int index = -1) override;
    void up() override;

    inline int level() const override
    {
        return static_cast<int>(nodes_.size()) - 1;
    }

    bool is_null_type() const override;
    bool get_bool_value() const override;
    std::int64_t get_int_value() const override;
    double get_double_value() const override;
    std::string get_string_value() const override;
    unsigned int get_node_type_index() const override;
    bool get_record_key_set(std::set<std::string> &) const override;
    bool is_effectively_null() const override;
    int get_array_size() const override;
    bool is_record_type() const override;
    bool is_array_type() const override;

    PathPython *clone() const override;

    std::string get_record_tag() const override;

    PathBase *find_ref_node() override;

protected:

    PathPython();

    PythonNodeType node_type() const;

    inline const py::object &head() const
    {
        return nodes_.back();
    }

    py::object root_node_;

    std::vector<py::object> nodes_;
};


/**
 * @brief Output operator for PathPython.
 *
 * Mainly for debugging purposes and error messages.
 */
std::ostream& operator<<(std::ostream& stream, const PathPython& path);


#pragma GCC visibility pop


} // namespace Input



#endif /* PATH_PYTHON_HH_ */
