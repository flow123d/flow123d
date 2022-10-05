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
 * @file    python_field_base.cc
 * @brief
 */

#include <pybind11.h>
#include <embed.h>    // everything needed for embedding
#include <stl.h>      // type conversion
#include "fields/python_field_base.hh"

PYBIND11_MODULE(flowpy, m) {
    m.doc() = "pybind11 Flow123D plugin"; // optional module docstring

    py::class_<PythonFieldBase>(m, "PythonFieldBaseCPP")
        .def(py::init<>())
        .def("_set_dict", &PythonFieldBase::set_dict)
        .def("_set_result", &PythonFieldBase::set_result)
        .def("_add_to_dict", &PythonFieldBase::add_to_dict)
//        .def("_set_result_data", &PythonFieldBase::set_result_data)
//        .def("_add_to_dict_data", &PythonFieldBase::add_to_dict_data)
        .def("_print_fields", &PythonFieldBase::print_fields)
        .def("_print_result", &PythonFieldBase::print_result)
        .def_property("t", &PythonFieldBase::get_time, &PythonFieldBase::set_time)
        .def_property("result", &PythonFieldBase::get_field_result, &PythonFieldBase::set_field_result)
        .def_property("f_dict", &PythonFieldBase::get_fields_dict, &PythonFieldBase::set_fields_dict);

    py::class_<FieldCacheProxy>(m, "FieldCacheProxy")
        .def(py::init<std::string, std::vector<ssize_t>, std::vector<double> >())
        .def("field_name", &FieldCacheProxy::field_name)
        .def("n_rows", &FieldCacheProxy::n_rows)
        .def("n_cols", &FieldCacheProxy::n_cols)
        .def("field_cache_ptr", &FieldCacheProxy::field_cache_ptr);

}
