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
 * @file    python_field_base.hh
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
    : field_name_(field_name), n_comp_(n_comp), field_cache_ptr_(field_cache_ptr)
    {}

    /// Copy constructor
    FieldCacheProxy(const FieldCacheProxy &other)
    : field_name_(other.field_name_), n_comp_(other.n_comp_), field_cache_ptr_(other.field_cache_ptr_)
    {}

    /// Getters
    const std::string &field_name() const { return field_name_; }
    const ssize_t &n_comp() const { return n_comp_; }
    std::vector<double> &field_cache_ptr() { return field_cache_ptr_; }
private:
    std::string field_name_;
    ssize_t n_comp_;
    std::vector<double> field_cache_ptr_;
};


class PythonFieldBase
{
protected:
    py::array create_array_with_data(std::vector<double> &data, ssize_t n_comp)
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

public:
    PythonFieldBase()
    {}

    PythonFieldBase(std::vector<FieldCacheProxy> &data, FieldCacheProxy &result)
    {
        py::dtype d_type("float64");
        set_dict(data, result);
    }

    void set_dict(std::vector<FieldCacheProxy> &data, FieldCacheProxy &result)
    {
        // Fill dictionary of input fields
    	fields_dict_.clear();
        for (uint i=0; i<data.size(); ++i) {
            fields_dict_[data[i].field_name().c_str()] =
                    create_array_with_data(data[i].field_cache_ptr(), data[i].n_comp());
        }
        // Fill array of result field
        {
            fields_dict_[result.field_name().c_str()] =
                    create_array_with_data(result.field_cache_ptr(), result.n_comp());
            field_result_ = result.field_name();
        }
    }

    void set_time(double t)
    {
        this->time_ = t;
    }

    double get_time() const
    {
        return this->time_;
        // set self.t in Python
    }

    py::dict &get_fields_dict()
    {
        return fields_dict_;
    }

    void set_fields_dict(py::dict &dict)
    {
        this->fields_dict_ = dict;
    }

    py::array get_field_result()
    {
        return this->fields_dict_[field_result_.c_str()];
    }

    void set_field_result(py::array &res)
    {
        this->fields_dict_[field_result_.c_str()] = res;
    }

//    void set_result_data(std::string field_name, double *data, ssize_t n_rows, ssize_t n_cols, ssize_t size)
//    {
//        fields_dict_[field_name.c_str()] = this->create_array_with_data(data, n_rows*n_cols, size);
//        field_result_ = field_name;
//    }
//
//    void add_to_dict_data(std::string field_name, std::vector<double> &data, ssize_t n_comp)
//    {
//        fields_dict_[field_name.c_str()] = this->create_array_with_data(data, n_comp);
//    }

    void print_fields() const
    {
        std::cout << "Dictionary contains fields: " << std::endl;
        for (auto &item : fields_dict_)
        {
            std::cout << " - " << item.first << ":" << std::endl << item.second << std::endl;
        }
    }

    void print_result() const
    {
        std::cout << "Result is: " << field_result_ << std::endl;
        for (auto item : fields_dict_)
        {
            std::string field_name = item.first.attr("__str__")().cast<std::string>();
            if (field_name == field_result_) std::cout << item.second << std::endl;
        }
    }

protected:
    py::dict fields_dict_;
    std::string field_result_;
    double time_;
};

#pragma GCC visibility pop

#endif /* PYTHON_FIELD_BASE_HH_ */
