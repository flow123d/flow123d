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
 * @file    field_python.impl.hh
 * @brief   
 */

#ifndef FIELD_PYTHON_IMPL_HH_
#define FIELD_PYTHON_IMPL_HH_


#include <type_traits>
#include "fields/field_python.hh"
#include "fields/field_set.hh"
#include "python/flowpy/python_field_base.hh"
#include <include/pybind11/pybind11.h>
#include <include/pybind11/eval.h>

namespace py = pybind11;

/// Implementation.

namespace it = Input::Type;

FLOW123D_FORCE_LINK_IN_CHILD(field_python)



template <int spacedim, class Value>
const Input::Type::Record & FieldPython<spacedim, Value>::get_input_type()
{
    return it::Record("FieldPython", FieldAlgorithmBase<spacedim,Value>::template_name()+" Field given by a Python script.")
		.derive_from(FieldAlgorithmBase<spacedim, Value>::get_input_type())
		.copy_keys(FieldAlgorithmBase<spacedim, Value>::get_field_algo_common_keys())
		.declare_key("source_file", it::String(), it::Default::obligatory(),
				"Python script given as external file in format 'dir'.'file_name' without .py extension")
		.declare_key("class", it::String(), it::Default::obligatory(),
				"Function in the given script that returns tuple containing components of the return type.\n"
				"For NxM tensor values: tensor(row,col) = tuple( M*row + col ).")
		//.declare_key("units", FieldAlgorithmBase<spacedim, Value>::get_field_algo_common_keys(), it::Default::optional(),
		//		"Definition of unit.")
		.close();
}


template <int spacedim, class Value>
const int FieldPython<spacedim, Value>::registrar =
		Input::register_class< FieldPython<spacedim, Value>, unsigned int >("FieldPython") +
		FieldPython<spacedim, Value>::get_input_type().size();



template <int spacedim, class Value>
FieldPython<spacedim, Value>::FieldPython(unsigned int n_comp)
: FieldAlgorithmBase<spacedim, Value>( n_comp)
{
	this->is_constant_in_space_ = false;
}



template <int spacedim, class Value>
void FieldPython<spacedim, Value>::set_python_field_from_string(const string &python_source, const string &func_name)
{
    p_module_ = PythonLoader::load_module_from_string("python_field_"+func_name, func_name, python_source);
    set_func("python_field_"+func_name, func_name);
}





template <int spacedim, class Value>
void FieldPython<spacedim, Value>::init_from_input(const Input::Record &rec, const struct FieldAlgoBaseInitData& init_data) {
	this->init_unit_conversion_coefficient(rec, init_data);
	this->field_name_ = init_data.field_name_;

    std::string source_file = rec.val<string>("source_file");
    std::string source_class = rec.val<string>("class");
    try {
        set_python_field_from_class( source_file, source_class );
    } INPUT_CATCH(FilePath::ExcFileOpen, FilePath::EI_Address_String, rec)

    in_rec_ = rec;
}



template <int spacedim, class Value>
void FieldPython<spacedim, Value>::set_python_field_from_file(const string &file_name, const string &func_name)
{
    p_module_ = PythonLoader::load_module_from_file( string(file_name) );
    set_func(file_name, func_name);
}




template <int spacedim, class Value>
void FieldPython<spacedim, Value>::set_python_field_from_class(const string &file_name, const string &class_name)
{
    internal::PythonWrapper::initialize();

    p_module_ = PythonLoader::load_module_from_file( string(file_name) );
    try {
        p_class_ = p_module_.attr(class_name.c_str());
    } catch (const py::error_already_set &ex) {
        PythonLoader::throw_error(ex);
    }
}




template <int spacedim, class Value>
void FieldPython<spacedim, Value>::set_func(const string &module_name, const string &func_name)
{
    try {
        p_func_ = p_module_.attr(func_name.c_str());
    } catch (const py::error_already_set &ex) {
        PythonLoader::throw_error(ex);
    }

    // try field call
    ASSERT_PERMANENT_EQ(spacedim, 3);
    double x=1.0, y=2.0, z=3.0;
    try {
        p_value_ = p_func_(x, y, z);
    } catch (const py::error_already_set &e) {
        PythonLoader::throw_error(e);
    }

    if ( ! PyTuple_Check( p_value_.ptr()) ) {
        stringstream ss;
        ss << "Field '" << func_name << "' from the python module: " << module_name << " doesn't return Tuple." << endl;
        THROW( ExcMessage() << EI_Message( ss.str() ));
    }

    unsigned int size = p_value_.size();

    unsigned int value_size=this->value_.n_rows() * this->value_.n_cols();
    if ( size !=  value_size) {
        THROW( ExcInvalidCompNumber() << EI_FuncName(func_name) << EI_PModule(module_name) << EI_Size(size) << EI_ValueSize(value_size) );
    }

}

/**
 * Returns one value in one given point. ResultType can be used to avoid some costly calculation if the result is trivial.
 */
template <int spacedim, class Value>
typename Value::return_type const & FieldPython<spacedim, Value>::value(const Point &p, const ElementAccessor<spacedim> &elm)
{
    set_value(p,elm, this->value_);
    this->value_.scale(this->unit_conversion_coefficient_);
    return this->r_value_;
}


/**
 * Returns std::vector of scalar values in several points at once.
 */
template <int spacedim, class Value>
void FieldPython<spacedim, Value>::value_list (const Armor::array &point_list, const ElementAccessor<spacedim> &elm,
                   std::vector<typename Value::return_type>  &value_list)
{
	ASSERT_EQ( point_list.size(), value_list.size() );
    ASSERT( point_list.n_rows() == spacedim && point_list.n_cols() == 1 ).error("Invalid point size.\n");
    for(unsigned int i=0; i< point_list.size(); i++) {
        Value envelope(value_list[i]);
        ASSERT_EQ( envelope.n_rows(), this->value_.n_rows() )(i)
                .error("value_list[i] has wrong number of rows\n");
        set_value(point_list.vec<spacedim>(i), elm, envelope );
        envelope.scale(this->unit_conversion_coefficient_);
    }
}

/**
* Returns one vector value in one given point.
*/
template <int spacedim, class Value>
void FieldPython<spacedim, Value>::set_value(FMT_UNUSED const Point &p, FMT_UNUSED const ElementAccessor<spacedim> &elm, FMT_UNUSED Value &value)
{
    p_value_ = p_func_(p[0], p[1], p[2]);

    unsigned int pos =0;
    for(unsigned int row=0; row < value.n_rows(); row++)
        for(unsigned int col=0; col < value.n_cols(); col++, pos++)
            if ( std::is_integral< typename Value::element_type >::value ) value(row,col) = p_value_[pos].cast<int>();
            else value(row,col) = p_value_[pos].cast<double>();
}




template <int spacedim, class Value>
std::vector<const FieldCommon * > FieldPython<spacedim, Value>::set_dependency(FieldSet &field_set) {
	std::vector<const FieldCommon *> required_fields;
    py::list field_list;

    try {
        p_func_ = p_class_.attr("used_fields");
        field_list = p_func_(p_class_);
    } catch (const py::error_already_set &ex) {
        PythonLoader::throw_error(ex);
    }

    std::vector<FieldCacheProxy> field_data;
    std::vector<ssize_t> field_shape(2);

    for (auto f : field_list) {
        std::string field_name = f.cast<std::string>();
        auto field_ptr = field_set.field(field_name);
        if (field_ptr != nullptr) required_fields.push_back( field_ptr );
        else {
            field_ptr = field_set.user_field(field_name, this->time_);
            if (field_ptr != nullptr) required_fields.push_back( field_ptr );
            else THROW( ExcUnknownField() << EI_Field(field_name) << Input::EI_Address( in_rec_.address_string() ) );
        }
        field_shape[0] = ssize_t(field_ptr->shape_[0]);
        field_shape[1] = ssize_t(field_ptr->shape_[1]);
        field_data.emplace_back(field_name, field_shape, field_ptr->value_cache()->data_);
    }

    auto self_ptr = field_set.field(this->field_name_); // instance of FieldCommon of this field
    std::vector<ssize_t> result_shape = { Value::NRows_, Value::NCols_ };
    double * cache_data = self_ptr->value_cache()->data_;
    std::vector<double> cache_vec(cache_data, cache_data+CacheMapElementNumber::get());
    FieldCacheProxy result_data(this->field_name_, result_shape, cache_vec);
    p_func_ = p_class_.attr("set_dict");
    py::list field_data_list = py::cast(field_data);
    p_func_(p_class_, field_data_list, result_data);
    return required_fields;
}



template <int spacedim, class Value>
void FieldPython<spacedim, Value>::cache_update(FMT_UNUSED FieldValueCache<typename Value::element_type> &data_cache,
        ElementCacheMap &cache_map, unsigned int region_patch_idx)
{
    unsigned int reg_chunk_begin = cache_map.region_chunk_begin(region_patch_idx);
    unsigned int reg_chunk_end = cache_map.region_chunk_end(region_patch_idx);
    p_func_ = p_class_.attr("cache_update");
    p_func_(p_class_, reg_chunk_begin, reg_chunk_end);
}



template <int spacedim, class Value>
FieldPython<spacedim, Value>::~FieldPython() {}


namespace internal {

PythonWrapper::PythonWrapper() {
    std::string fname = std::string(FLOW123D_SOURCE_DIR) + "/src/python/flowpy/python_field_base.py";
    FilePath flowpy_path(fname, FilePath::input_file);
    std::string parent_path = flowpy_path.parent_path(); // add path to PythonPath
    PythonLoader::add_sys_path(parent_path);
    py::module_ flowpy_module = py::module_::import("flowpy");
    py::module_ pyfield_module = py::module_::import("python_field_base");
    flowpy_module.add_object("PythonFieldBase", pyfield_module.attr("PythonFieldBase"));
}

PythonWrapper &PythonWrapper::initialize() {
    static PythonWrapper python_wrapper;
    return python_wrapper;
}

}


#endif /* FIELD_PYTHON_IMPL_HH_ */
