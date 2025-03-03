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

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/numpy.h>
//#include <pybind11/common.h>
#include "fem/element_cache_map.hh"

namespace py = pybind11;
#pragma GCC visibility push(hidden)

/// Helper class, holds data of one field
class FieldCacheProxy
{
public:
	/**
	 * Method encapsulates FieldValueCache data array for usage in Python.
	 * Allows to create C++ and Python objects above shared block of memory.
	 */
	static py::buffer_info field_proxy_get_buffer(FieldCacheProxy &proxy)
	{
	    ssize_t              n_comp  = ( (proxy.shape_.size()==1) ? proxy.shape_[0] : proxy.shape_[0]*proxy.shape_[1]);
	    ssize_t              size    = proxy.data_size_ / n_comp;
	    std::vector<ssize_t> shape;
	    std::vector<ssize_t> strides;

	    if (proxy.shape_[0] > 1) { // add dimensions only for vector and tensor
	        shape.push_back(proxy.shape_[0]);
	        if (proxy.shape_.size() == 2) shape.push_back(proxy.shape_[1]);
	    }
	    shape.push_back(size);

	    ssize_t n_dim = shape.size();
	    strides.resize(n_dim);
	    strides[n_dim-1] = sizeof(double);
	    for(uint i=n_dim-1; i>0; i--) {
	        strides[i-1] = strides[i] * shape[i];
	    }

	    // create n_dim NumPy array
	    return  py::buffer_info(
	    	proxy.field_cache_data_,                 /* data as contiguous array  */
	        sizeof(double),                          /* size of one scalar        */
	        py::format_descriptor<double>::format(), /* data type                 */
			n_dim,                                   /* number of dimensions      */
	        shape,                                   /* shape of the matrix       */
	        strides                                  /* strides for each axis     */
	    );
	}

    /// Constructor
    FieldCacheProxy(std::string field_name, std::vector<uint> shape, double * field_cache_data, uint data_size)
    : field_name_(field_name), shape_(shape), field_cache_data_(field_cache_data), data_size_(data_size)
    {}

    /// Copy constructor
    FieldCacheProxy(const FieldCacheProxy &other)
    : field_name_(other.field_name_), shape_(other.shape_), field_cache_data_(other.field_cache_data_), data_size_(other.data_size_)
    {}

    /// Getter returns field name
    const std::string &field_name() const { return field_name_; }


private:
    std::string field_name_;
    std::vector<uint> shape_;
    double *field_cache_data_;
    uint data_size_;
};

#pragma GCC visibility pop

#endif /* PYTHON_FIELD_BASE_HH_ */
