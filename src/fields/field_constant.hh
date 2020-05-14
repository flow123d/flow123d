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
 * @file    field_constant.hh
 * @brief   
 */

#ifndef FIELD_CONSTANT_HH_
#define FIELD_CONSTANT_HH_

#include <stdio.h>                            // for sprintf
#include <string.h>                           // for memcpy
#include <boost/exception/info.hpp>           // for operator<<, error_info:...
#include <new>                                // for operator new[]
#include <string>                             // for operator<<
#include <vector>                             // for vector
#include <armadillo>
#include "fields/field_algo_base.hh"          // for FieldAlgorithmBase
#include "fields/field_values.hh"             // for FieldValue<>::TensorFixed
#include "tools/unit_converter.hh"           // for ExcInvalidUnit, ExcNonc...
#include "input/accessors.hh"                 // for ExcTypeMismatch, ExcAcc...
#include "input/accessors_impl.hh"            // for Record::val, Iterator::...
#include "input/storage.hh"                   // for ExcStorageTypeMismatch
#include "input/type_record.hh"               // for Record::ExcRecordKeyNot...
#include "system/exceptions.hh"               // for ExcAssertMsg::~ExcAsser...
template <int spacedim> class ElementAccessor;


/**
 * Class representing spatially constant fields.
 *
 */
template <int spacedim, class Value>
class FieldConstant : public FieldAlgorithmBase<spacedim, Value>
{
public:
    typedef typename FieldAlgorithmBase<spacedim, Value>::Point Point;
    typedef FieldAlgorithmBase<spacedim, Value> FactoryBaseType;

    /**
     * Default constructor, optionally we need number of components @p n_comp in the case of Vector valued fields.
     */
    FieldConstant(unsigned int n_comp=0);


    /**
     * Return Record for initialization of FieldConstant that is derived from Abstract given by @p a_type
     * and the individual elements of the possible Value (vector, tensor) have Input::Type @p eit.
     */
    static const Input::Type::Record & get_input_type();

    /**
     * Smart setter from the given value to return.
     */
    FieldConstant<spacedim, Value> &set_value(const typename Value::return_type &val);

    /**
     * This method initialize actual value of the field given from the given Input::Record @p rec.
     */
    virtual void init_from_input(const Input::Record &rec, const struct FieldAlgoBaseInitData& init_data) override;



    /**
     * Returns one value in one given point. ResultType can be used to avoid some costly calculation if the result is trivial.
     */
    virtual typename Value::return_type const &value(const Point &p, const ElementAccessor<spacedim> &elm) override;

    /**
     * Returns std::vector of scalar values in several points at once.
     */
    virtual void value_list (const Armor::array &point_list, const ElementAccessor<spacedim> &elm,
                       std::vector<typename Value::return_type>  &value_list) override;

    void cache_update(FieldValueCache<typename Value::element_type> &data_cache,
			ElementCacheMap &cache_map, unsigned int region_idx) override;


    virtual ~FieldConstant();

protected:
    /// Compare field value with given minimal and maximal limits.
    void check_field_limits(const Input::Record &rec, const struct FieldAlgoBaseInitData& init_data);

private:
    /// Registrar of class to factory
    static const int registrar;

};


#endif /* FIELD_CONSTANT_HH_ */
