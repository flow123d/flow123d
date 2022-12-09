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
#include <pybind11.h>

using namespace std;
namespace py = pybind11;

// Pybind11 needs set visibility to hidden (see https://pybind11.readthedocs.io/en/stable/faq.html).
#pragma GCC visibility push(hidden)

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

    TYPEDEF_ERR_INFO( EI_FuncName, std::string);
    TYPEDEF_ERR_INFO( EI_PModule, std::string);
    TYPEDEF_ERR_INFO( EI_Size, unsigned int);
    TYPEDEF_ERR_INFO( EI_ValueSize, unsigned int);
    DECLARE_EXCEPTION( ExcInvalidCompNumber, << "Field " << EI_FuncName::qval << " from the python module: " << EI_PModule::val
            << " returns " << EI_Size::val << " components but should return " << EI_ValueSize::val << " components.\n" );


    FieldPython(unsigned int n_comp=0);

    virtual void init_from_input(const Input::Record &rec, const struct FieldAlgoBaseInitData& init_data);

    static const Input::Type::Record & get_input_type();

    /**
     * Set the source in a string and name of the field to be called.
     */
    void set_python_field_from_class(const string &file_name, const string &class_name);

    /**
     * Overload @p FieldAlgorithmBase::cache_reinit
     *
     * Reinit dictionary of used fields and update result field to dictionary of resul fields.
     */
    void cache_reinit(const ElementCacheMap &cache_map) override;

    void cache_update(FieldValueCache<typename Value::element_type> &data_cache,
			ElementCacheMap &cache_map, unsigned int region_patch_idx) override;

    /**
     * Returns list of fields on which this field depends.
     */
    std::vector<const FieldCommon *> set_dependency(FieldSet &field_set) override;

    virtual ~FieldPython();

private:
    /// Registrar of class to factory
    static const int registrar;

    /// Accessor to Input::Record
    Input::Record in_rec_;

    /// Field name is necessary for set result
    std::string field_name_;

    /// Holds python instance of user class.
    py::object user_class_instance_;

    /// List of fields on which this field depends
	std::vector<const FieldCommon * > required_fields_;

    /// Pointer to FieldCommon that holds this fields (stores in set_dependency and uses in cache_reinit)
	const FieldCommon * self_field_ptr_;

};


#pragma GCC visibility pop

#endif /* FUNCTION_PYTHON_HH_ */
