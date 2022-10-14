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

    py::class_<FieldCacheProxy>(m, "FieldCacheProxy")
        .def(py::init<std::string, ssize_t, std::vector<double> >())
        //.def_property("field_name", &FieldCacheProxy::get_field_name, &FieldCacheProxy::set_field_name)
        .def("field_name", &FieldCacheProxy::field_name)
        .def("field_cache_array", &FieldCacheProxy::field_cache_array);

}
