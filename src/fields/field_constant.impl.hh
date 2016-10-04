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
 * @file    field_constant.impl.hh
 * @brief   
 */

#ifndef FIELD_CONSTANT_IMPL_HH_
#define FIELD_CONSTANT_IMPL_HH_

#include "fields/field_constant.hh"
#include "input/input_type.hh"


/// Implementation.

namespace it = Input::Type;

FLOW123D_FORCE_LINK_IN_CHILD(field_constant)


template <int spacedim, class Value>
const Input::Type::Record & FieldConstant<spacedim, Value>::get_input_type()
{
    return it::Record("FieldConstant", FieldAlgorithmBase<spacedim,Value>::template_name()+" Field constant in space.")
        .derive_from(FieldAlgorithmBase<spacedim, Value>::get_input_type())
        .declare_key("value", Value::get_input_type(), it::Default::obligatory(),
                                    "Value of the constant field.\n"
                                    "For vector values, you can use scalar value to enter constant vector.\n"
                                    "For square (($N\\times N$))-matrix values, you can use:\n"
                                    " - vector of size (($N$)) to enter diagonal matrix\n\n"
                                    " - vector of size (($\\frac12N(N+1)$)) to enter symmetric matrix (upper triangle, row by row)\n"
                                    " - scalar to enter multiple of the unit matrix." )
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
{}


template <int spacedim, class Value>
FieldConstant<spacedim, Value> &FieldConstant<spacedim, Value>::set_value(const typename Value::return_type &val)
{
    this->r_value_ = val;

    return *this;
}


template <int spacedim, class Value>
void FieldConstant<spacedim, Value>::init_from_input(const Input::Record &rec) {
    this->value_.init_from_input( rec.val<typename Value::AccessType>("value") );

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



/**
 * Returns one value in one given point. ResultType can be used to avoid some costly calculation if the result is trivial.
 */
template <int spacedim, class Value>
typename Value::return_type const & FieldConstant<spacedim, Value>::value(const Point &p, const ElementAccessor<spacedim> &elm)
{
    return this->r_value_;
}



/**
 * Returns std::vector of scalar values in several points at once.
 */
template <int spacedim, class Value>
void FieldConstant<spacedim, Value>::value_list (const std::vector< Point >  &point_list, const ElementAccessor<spacedim> &elm,
                   std::vector<typename Value::return_type>  &value_list)
{
	OLD_ASSERT_EQUAL( point_list.size(), value_list.size() );

    for(unsigned int i=0; i< point_list.size(); i++) {
    	OLD_ASSERT( Value(value_list[i]).n_rows()==this->value_.n_rows(),
                "value_list[%d] has wrong number of rows: %d; should match number of components: %d\n",
                i, Value(value_list[i]).n_rows(),this->value_.n_rows());


        value_list[i]=this->r_value_;
    }
}



template <int spacedim, class Value>
FieldConstant<spacedim, Value>::~FieldConstant() {
}



#endif /* FIELD_CONSTANT_IMPL_HH_ */
