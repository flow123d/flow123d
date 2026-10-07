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
 * @file    path_python.cc
 * @brief
 */

#include "input/path_python.hh"
#include "input/reader_internal_base.hh"
#include "input/reader_to_storage.hh"
#include "system/system.hh"


namespace Input {
using namespace std;


PathPython::PathPython(const py::dict &root)
: PathPython()
{
    root_node_ = root;
    nodes_.push_back(root_node_);
}

PathPython::PathPython()
: PathBase()
{
    json_type_names.push_back("Python dict");
    json_type_names.push_back("Python list");
    json_type_names.push_back("Python string");
    json_type_names.push_back("Python bool");
    json_type_names.push_back("Python int");
    json_type_names.push_back("Python real");
    json_type_names.push_back("Python None");
    json_type_names.push_back("other Python type");
    json_type_names.push_back("undefined type");
}

PathPython::~PathPython() = default;

bool PathPython::down(unsigned int index) {
    const py::object &head_node = nodes_.back();

    if (!is_array_type()) {
        return false;
    }

    py::list array = head_node.cast<py::list>();

    if (index >= array.size()) {
        return false;
    }

    path_.push_back(
        std::make_pair(static_cast<int>(index), std::string(""))
    );

    nodes_.push_back(array[index]);

    return true;
}

bool PathPython::down(const std::string &key, int index) {
    const py::object &head_node = nodes_.back();

    if (!is_record_type()) {
        return false;
    }

    py::dict dict = head_node.cast<py::dict>();

    py::str py_key(key);

    if (!dict.contains(py_key)) {
        return false;
    }

    path_.push_back(std::make_pair(index, key));
    nodes_.push_back(dict[py_key]);

    return true;
}

void PathPython::up() {
    if (path_.size() > 1) {
        path_.pop_back();
        nodes_.pop_back();
    }
}

unsigned int PathPython::get_node_type_index() const {
    const py::object &node = head();

    if (node.is_none()) {
        return ValueTypes::null_type;
    }

    if (py::isinstance<py::dict>(node)) {
        return ValueTypes::obj_type;
    }

    if (py::isinstance<py::list>(node)) {
        return ValueTypes::array_type;
    }

    if (py::isinstance<py::str>(node)) {
        return ValueTypes::str_type;
    }

    if (py::isinstance<py::bool_>(node)) {
        return ValueTypes::bool_type;
    }

    if (py::isinstance<py::int_>(node)) {
        return ValueTypes::int_type;
    }

    if (py::isinstance<py::float_>(node)) {
        return ValueTypes::real_type;
    }

    return ValueTypes::scalar_type;
}

bool PathPython::is_record_type() const {
    return py::isinstance<py::dict>(head());
}

bool PathPython::is_array_type() const {
    return py::isinstance<py::list>(head());
}

bool PathPython::get_bool_value() const {
    if (py::isinstance<py::bool_>(head())) {
        return head().cast<bool>();
    } else {
        THROW(ReaderInternalBase::ExcInputError());
    }

    return false;
}

std::int64_t PathPython::get_int_value() const {
    if (py::isinstance<py::int_>(head())) {
        return head().cast<std::int64_t>();
    } else {
        THROW(ReaderInternalBase::ExcInputError());
    }

    return 0;
}

double PathPython::get_double_value() const {
    if (py::isinstance<py::float_>(head()) ||
        py::isinstance<py::int_>(head())) {

        return head().cast<double>();
    } else {
        THROW(ReaderInternalBase::ExcInputError());
    }

    return 0.0;
}

std::string PathPython::get_string_value() const {
    if (py::isinstance<py::str>(head())) {
        return head().cast<std::string>();
    } else {
        THROW(ReaderInternalBase::ExcInputError());
    }

    return "";
}

int PathPython::get_array_size() const {
    if (!is_array_type()) {
        return -1;
    }

    return static_cast<int>(
        head().cast<py::list>().size()
    );
}

bool PathPython::is_null_type() const {
    return head().is_none();
}

bool PathPython::get_record_key_set(std::set<std::string> &keys_list) const {
    if (!is_record_type()) {
        return false;
    }

    py::dict dict = head().cast<py::dict>();

    for (auto item : dict) {
        py::handle key = item.first;

        if (!py::isinstance<py::str>(key)) {
            THROW(ReaderInternalBase::ExcInputError());
        }

        keys_list.insert(key.cast<std::string>());
    }

    return true;
}

bool PathPython::is_effectively_null() const {
    return false;
}

std::string PathPython::get_record_tag() const {
//    // Variant with support of REF
//    std::string tag_value;
//
//    if (is_record_type()) {
//        PathPython type_path(*this);
//
//        if (type_path.down("TYPE")) {
//            PathBase *ref_path = type_path.find_ref_node();
//
//            if (ref_path) {
//                tag_value = ref_path->get_string_value();
//                delete ref_path;
//            } else {
//                tag_value = type_path.get_string_value();
//            }
//        }
//    }
//
//    return tag_value;

    // Variant without support of REF
    if (!is_record_type()) {
        return "";
    }

    PathPython type_path(*this);

    if (!type_path.down("TYPE")) {
        return "";
    }

    return type_path.get_string_value();
}

PathBase * PathPython::find_ref_node() {
    // REF is not supported yet.
    return nullptr;
}

PathPython *PathPython::clone() const {
    return new PathPython(*this);
}

std::ostream& operator<<(std::ostream& stream, const PathPython& path) {
    path.output(stream);
    return stream;
}


} // namespace Input
