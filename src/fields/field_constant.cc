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
 * @file    field_constant.cc
 * @brief   
 */

#include "fields/field_constant.hh"
#include "fields/field_algo_base.impl.hh"
#include "fields/field_instances.hh"	// for instantiation macros
#include "input/input_type.hh"
#include "system/armor.hh"


/// Implementation.

namespace it = Input::Type;

FLOW123D_FORCE_LINK_IN_CHILD(field_constant)


template <int spacedim, class Value>
const Input::Type::Record & FieldConstant<spacedim, Value>::get_input_type()
{
    return it::Record("FieldConstant", FieldAlgorithmBase<spacedim,Value>::template_name()+" Field constant in space.")
        .derive_from(FieldAlgorithmBase<spacedim, Value>::get_input_type())
        .copy_keys(FieldAlgorithmBase<spacedim, Value>::get_field_algo_common_keys())
        .declare_key("value", Value::get_input_type(), it::Default::obligatory(),
                                    "Value of the constant field. "
                                    "For vector values, you can use scalar value to enter constant vector. "
                                    "For square (($N\\times N$))-matrix values, you can use: "
                                    " - vector of size (($N$)) to enter diagonal matrix\n\n"
                                    " - vector of size (($\\frac12N(N+1)$)) to enter symmetric matrix (upper triangle, row by row)\n"
                                    " - scalar to enter multiple of the unit matrix." )
		//.declare_key("unit", FieldAlgorithmBase<spacedim, Value>::get_input_type_unit_si(), it::Default::optional(),
		//							"Definition of unit.")
        .allow_auto_conversion("value")
		.close();
}

template <int spacedim, class Value>
const int FieldConstant<spacedim, Value>::registrar =
		Input::register_class< FieldConstant<spacedim, Value>, unsigned int >("FieldConstant") +
		FieldConstant<spacedim, Value>::get_input_type().size();


template <int spacedim, class Value>
FieldConstant<spacedim, Value>::FieldConstant( unsigned int n_comp)
: FieldAlgorithmBase<spacedim, Value>(n_comp)
{
	this->is_constant_in_space_ = true;
}


template <int spacedim, class Value>
FieldConstant<spacedim, Value> &FieldConstant<spacedim, Value>::set_value(const typename Value::return_type &val)
{
    this->r_value_ = val;

    return *this;
}


template <int spacedim, class Value>
void FieldConstant<spacedim, Value>::init_from_input(const Input::Record &rec, const struct FieldAlgoBaseInitData& init_data) {
	this->init_unit_conversion_coefficient(rec, init_data);


    this->value_.init_from_input( rec.val<typename Value::AccessType>("value") );
    this->value_.scale(this->unit_conversion_coefficient_);
    this->check_field_limits(rec, init_data);

    typename Value::return_type tmp_value;
    Value tmp_field_value(tmp_value);
    tmp_field_value.set_n_comp(this->n_comp());

    tmp_field_value.zeros();
    if ( this->value_.equal_to(tmp_value) ) {
        this->field_result_ = result_zeros;
        return;
    }


    tmp_field_value.ones();
    if ( this->value_.equal_to(tmp_value) ) {
        this->field_result_ = result_ones;
        return;
    }

    // This check must be the last one, since for scalar and vector values ones() == eye().
    // For vector, eye() does nothing. So, the value of tmp_value remains equal to ones().
    tmp_field_value.eye();
    if ( this->value_.equal_to(tmp_value) ) {
        this->field_result_ = result_eye;
        return;
    }


    this->field_result_ = result_constant;
}



template <int spacedim, class Value>
void FieldConstant<spacedim, Value>::cache_update(FieldValueCache<typename Value::element_type> &data_cache,
		ElementCacheMap &cache_map, unsigned int region_patch_idx)
{
    unsigned int reg_chunk_begin = cache_map.region_chunk_begin(region_patch_idx);
    unsigned int reg_chunk_end = cache_map.region_chunk_end(region_patch_idx);
    Armor::ArmaMat<typename Value::element_type, Value::NRows_, Value::NCols_> mat_value( const_cast<typename Value::element_type*>(this->value_.mem_ptr()) );
    for (unsigned int i_cache = reg_chunk_begin; i_cache < reg_chunk_end; ++i_cache)
        data_cache.set(i_cache) = mat_value;
}


template <int spacedim, class Value>
void FieldConstant<spacedim, Value>::check_field_limits(const Input::Record &rec, const struct FieldAlgoBaseInitData& init_data)
{
    if (Value::is_scalable())
        for( unsigned int row=0; row<this->value_.n_rows(); row++)
            for( unsigned int col=0; col<this->value_.n_cols(); col++) {
            	if ( (this->value_(row,col) < init_data.limits_.first) || (this->value_(row,col) > init_data.limits_.second) ) {
                    WarningOut().fmt("Value '{}' of Field '{}' at address '{}' is out of limits: <{}, {}>\nUnit of the Field: [{}]\n",
                    		this->value_(row,col), init_data.field_name_, rec.address_string(),
							init_data.limits_.first, init_data.limits_.second, init_data.unit_si_.format_text() );
            	}
            }
}



template <int spacedim, class Value>
FieldConstant<spacedim, Value>::~FieldConstant() {
}


// Instantiations of FieldConstant
INSTANCE_ALL(FieldConstant)

// temporary solution for computing more fields at once in python
template class FieldConstant<3, FieldValue<0>::Vector >; // Necessary due to default value of the abstract.

