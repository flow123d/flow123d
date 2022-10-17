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
 * @file    python_field_proxy.hh
 * @brief
 */

#ifndef PYTHON_FIELD_BASE_HH_
#define PYTHON_FIELD_BASE_HH_

#include <pybind11.h>
#include <stl.h>
#include <numpy.h>
#include <detail/common.h>
#include "fields/field_value_cache.hh"

namespace py = pybind11;
#pragma GCC visibility push(hidden)

/// Helper class, holds data of one field
class FieldCacheProxy
{
public:
    /// Constructor
    FieldCacheProxy(std::string field_name, ssize_t n_comp, std::vector<double> field_cache_ptr)
    : field_name_(field_name)
    {
        field_cache_array_ = this->create_array_from_vec(field_cache_ptr, n_comp);
    }

    /// Copy constructor
    FieldCacheProxy(const FieldCacheProxy &other)
    : field_name_(other.field_name_), field_cache_array_(other.field_cache_array_)
    {}

    /// Getters
    const std::string &field_name() const { return field_name_; }
    py::array &field_cache_array() { return field_cache_array_; }
private:
    py::array create_array_from_vec(std::vector<double> &data, ssize_t n_comp)
    {
        ssize_t              size    = data.size() / n_comp;
        ssize_t              ndim    = 2;
        std::vector<ssize_t> shape   = { n_comp , size };
        std::vector<ssize_t> strides = { (long int)(sizeof(double)*size) , sizeof(double) };

        // create 2-D NumPy array
        return  py::array(py::buffer_info(
            &data[0],                                /* data as contiguous array  */
            sizeof(double),                          /* size of one scalar        */
            py::format_descriptor<double>::format(), /* data type                 */
            ndim,                                    /* number of dimensions      */
            shape,                                   /* shape of the matrix       */
            strides                                  /* strides for each axis     */
        ));
    }

    std::string field_name_;
    py::array field_cache_array_;
};


#pragma GCC visibility pop

#endif /* PYTHON_FIELD_BASE_HH_ */
