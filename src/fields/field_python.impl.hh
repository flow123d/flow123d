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

/// Implementation.

namespace it = Input::Type;

FLOW123D_FORCE_LINK_IN_CHILD(field_python)



template <int spacedim, class Value>
const Input::Type::Record & FieldPython<spacedim, Value>::get_input_type()
{
    return it::Record("FieldPython", FieldAlgorithmBase<spacedim,Value>::template_name()+" Field given by a Python script.")
		.derive_from(FieldAlgorithmBase<spacedim, Value>::get_input_type())
		.copy_keys(FieldAlgorithmBase<spacedim, Value>::get_field_algo_common_keys())
		.declare_key("script_string", it::String(), it::Default::read_time("Obligatory if 'script_file' is not given. "),
				"Python script given as in place string")
		.declare_key("script_file", it::FileName::input(), it::Default::read_time("Obligatory if 'script_striong' is not given. "),
				"Python script given as external file")
		.declare_key("function", it::String(), it::Default::obligatory(),
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

#ifdef FLOW123D_HAVE_PYTHON
    p_func_=NULL;
    p_module_=NULL;
    p_args_=NULL;
    p_value_=NULL;
#else
    THROW( ExcNoPythonSupport() );
#endif // FLOW123D_HAVE_PYTHON
}



template <int spacedim, class Value>
void FieldPython<spacedim, Value>::set_python_field_from_string(FMT_UNUSED const string &python_source, FMT_UNUSED const string &func_name)
{
#ifdef FLOW123D_HAVE_PYTHON
    p_module_ = PythonLoader::load_module_from_string("python_field_"+func_name, python_source);
    set_func(func_name);
#endif // FLOW123D_HAVE_PYTHON
}





template <int spacedim, class Value>
void FieldPython<spacedim, Value>::init_from_input(const Input::Record &rec, const struct FieldAlgoBaseInitData& init_data) {
	this->init_unit_conversion_coefficient(rec, init_data);

    Input::Iterator<string> it = rec.find<string>("script_string");
    if (it) {
        set_python_field_from_string( *it, rec.val<string>("function") );
    } else {
        Input::Iterator<FilePath> it = rec.find<FilePath>("script_file");
        if (! it) THROW( ExcNoPythonInit() );
        try {
            set_python_field_from_file( *it, rec.val<string>("function") );
        } INPUT_CATCH(FilePath::ExcFileOpen, FilePath::EI_Address_String, rec)
    }
}



template <int spacedim, class Value>
void FieldPython<spacedim, Value>::set_python_field_from_file(FMT_UNUSED const FilePath &file_name, FMT_UNUSED const string &func_name)
{
#ifdef FLOW123D_HAVE_PYTHON
    p_module_ = PythonLoader::load_module_from_file( string(file_name) );
    set_func(func_name);
#endif // FLOW123D_HAVE_PYTHON
}




template <int spacedim, class Value>
void FieldPython<spacedim, Value>::set_func(FMT_UNUSED const string &func_name)
{
#ifdef FLOW123D_HAVE_PYTHON
	p_func_ = PythonLoader::get_callable(p_module_, func_name);

    p_args_ = PyTuple_New( spacedim );

    // try field call
    for(unsigned int i = 0; i < spacedim; i++) {
        p_value_ = PyFloat_FromDouble( double(i) );
        PyTuple_SetItem(p_args_, i, p_value_);
    }
    p_value_ = PyObject_CallObject(p_func_, p_args_);
    PythonLoader::check_error();

    if ( ! PyTuple_Check( p_value_)) {
    	stringstream ss;
    	ss << "Field '" << func_name << "' from the python module: " << PyModule_GetName(p_module_) << " doesn't return Tuple." << endl;
        THROW( ExcMessage() << EI_Message( ss.str() ));
    }

    unsigned int size = PyTuple_Size( p_value_);

    unsigned int value_size=this->value_.n_rows() * this->value_.n_cols();
    if ( size !=  value_size) {
        THROW( ExcInvalidCompNumber() << EI_FuncName(func_name) << EI_PModule(PyModule_GetName(p_module_)) << EI_Size(size) << EI_ValueSize(value_size) );
    }

#endif // FLOW123D_HAVE_PYTHON

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
	OLD_ASSERT_EQUAL( point_list.size(), value_list.size() );
    ASSERT_DBG( point_list.n_rows() == spacedim && point_list.n_cols() == 1 ).error("Invalid point size.\n");
    for(unsigned int i=0; i< point_list.size(); i++) {
        Value envelope(value_list[i]);
        OLD_ASSERT( envelope.n_rows()==this->value_.n_rows(),
                "value_list[%d] has wrong number of rows: %d; should match number of components: %d\n",
                i, envelope.n_rows(),this->value_.n_rows());
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
#ifdef FLOW123D_HAVE_PYTHON
    for(unsigned int i = 0; i < spacedim; i++) {
        p_value_ = PyFloat_FromDouble( p[i] );
        PyTuple_SetItem(p_args_, i, p_value_);
    }
    p_value_ = PyObject_CallObject(p_func_, p_args_);
    PythonLoader::check_error();

    unsigned int pos =0;
    for(unsigned int row=0; row < value.n_rows(); row++)
        for(unsigned int col=0; col < value.n_cols(); col++, pos++)
            if ( std::is_integral< typename Value::element_type >::value ) value(row,col) = PyLong_AsLong( PyTuple_GetItem( p_value_, pos ) );
            else value(row,col) = PyFloat_AsDouble( PyTuple_GetItem( p_value_, pos ) );

#endif // FLOW123D_HAVE_PYTHON
}




template <int spacedim, class Value>
FieldPython<spacedim, Value>::~FieldPython() {
#ifdef FLOW123D_HAVE_PYTHON
    Py_CLEAR(p_module_);
    Py_CLEAR(p_func_);
    Py_CLEAR(p_value_);
    Py_CLEAR(p_args_);
#endif // FLOW123D_HAVE_PYTHON
}


#endif /* FIELD_PYTHON_IMPL_HH_ */
