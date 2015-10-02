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
 * @file    field_python.hh
 * @brief   
 * @todo    make FieldPython dummy class if we do not have python, so that we
 *          need not optional code elsewhere
 */

#ifndef FIELD_PYTHON_HH_
#define FIELD_PYTHON_HH_


#include "system/system.hh"
#include "system/python_loader.hh"
#include "fields/field_algo_base.hh"
#include "mesh/point.hh"
#include "input/factory.hh"

#include <string>
using namespace std;

/**
 *
 * This class assumes field python field with @p spacedim arguments containing coordinates of the given point.
 * The field should return  a tuple representing a vector value (possibly of size one for scalar fields)
 *
 * TODO:
 * - use rather only one argument - tuple representing the whole point
 * - time fields
 * - set some parameter in the python module
 * - for fields with one component allow python fields returning the value directly
 *
 */
template <int spacedim, class Value>
class FieldPython : public FieldAlgorithmBase<spacedim, Value>
{
public:
    typedef typename FieldAlgorithmBase<spacedim, Value>::Point Point;
    typedef FieldAlgorithmBase<spacedim, Value> FactoryBaseType;

    FieldPython(unsigned int n_comp=0);

    virtual void init_from_input(const Input::Record &rec);

    static const Input::Type::Record & get_input_type(Input::Type::AbstractRecord &a_type, const typename Value::ElementInputType *eit);

    /**
     * Set the file and field to be called.
     * TODO: use FilePath
     */
    void set_python_field_from_file( const FilePath &file_name, const string &func_name);

    /**
     * Set the source in a string and name of the field to be called.
     */
    void set_python_field_from_string( const string &python_source, const string &func_name);

    /**
     * Returns one value in one given point. ResultType can be used to avoid some costly calculation if the result is trivial.
     */
    virtual typename Value::return_type const &value(const Point &p, const ElementAccessor<spacedim> &elm);

    /**
     * Returns std::vector of scalar values in several points at once.
     */
    virtual void value_list (const std::vector< Point >  &point_list, const ElementAccessor<spacedim> &elm,
                       std::vector<typename Value::return_type>  &value_list);


    virtual ~FieldPython();

private:
    /// Registrar of class to factory
    static const int registrar;

    /**
     * Common part of set_python_field_from_* methods
     */
    void set_func(const string &func_name);

    /**
     * Implementation.
     */
    inline void set_value(const Point &p, const ElementAccessor<spacedim> &elm, Value &value);

#ifdef FLOW123D_HAVE_PYTHON
    PyObject *p_func_;
    PyObject *p_module_;
    mutable PyObject *p_args_;
    mutable PyObject *p_value_;
#endif // FLOW123D_HAVE_PYTHON

};



#endif /* FUNCTION_PYTHON_HH_ */
