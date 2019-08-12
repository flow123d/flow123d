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
 * @file    field_time_function.cc
 * @brief
 */

#include "fields/field_time_function.hh"
#include "fields/field_instances.hh"	// for instantiation macros
#include "fields/table_function.hh"

namespace it = Input::Type;

FLOW123D_FORCE_LINK_IN_CHILD(field_time_function)


template <int spacedim, class Value>
const Input::Type::Record & FieldTimeFunction<spacedim, Value>::get_input_type()
{
    return it::Record("FieldTimeFunction", FieldAlgorithmBase<spacedim,Value>::template_name()+" Field time-dependent function in space.")
        .derive_from(FieldAlgorithmBase<spacedim, Value>::get_input_type())
        .copy_keys(FieldAlgorithmBase<spacedim, Value>::get_field_algo_common_keys())
		.declare_key("time_function", TableFunction<Value>::get_input_type(), it::Default::obligatory(),
									"Values of time series initialization of Field.")
		//.declare_key("unit", FieldAlgorithmBase<spacedim, Value>::get_field_algo_common_keys(), it::Default::optional(),
		//							"Definition of unit.")
		//.declare_key("time_unit", it::String(), it::Default::read_time("Common unit of TimeGovernor."),
		//							"Definition of unit of all times defined in FieldTimeFunction.")
        .allow_auto_conversion("time_function")
		.close();
}

template <int spacedim, class Value>
const int FieldTimeFunction<spacedim, Value>::registrar =
		Input::register_class< FieldTimeFunction<spacedim, Value>, unsigned int >("FieldTimeFunction") +
		FieldTimeFunction<spacedim, Value>::get_input_type().size();


template <int spacedim, class Value>
FieldTimeFunction<spacedim, Value>::FieldTimeFunction( unsigned int n_comp)
: FieldConstant<spacedim, Value>(n_comp), unit_si_( UnitSI::dimensionless() )
{}


template <int spacedim, class Value>
void FieldTimeFunction<spacedim, Value>::init_from_input(const Input::Record &rec, const struct FieldAlgoBaseInitData& init_data)
{
	this->init_unit_conversion_coefficient(rec, init_data);
	this->field_name_ = init_data.field_name_;
	this->in_rec_ = rec;
	this->unit_si_ = init_data.unit_si_;
	this->limits_ = init_data.limits_;
}


template <int spacedim, class Value>
bool FieldTimeFunction<spacedim, Value>::set_time(const TimeStep &time)
{
	// Possible optimization: If we need set value_ repeatedly, we introduce TableFunction as class member and retrive it only once.
	TableFunction<Value> table_function;

	if (!Value::is_scalable()) {
		WarningOut().fmt("Setting key 'time_function' of non-floating point field at address {}\nValues will be skipped.\n",
				in_rec_.address_string());
	}
	table_function.init_from_input( in_rec_.val<Input::Record>("time_function"), time );
	this->r_value_ = table_function.value( time.end() );
    this->value_.scale(this->unit_conversion_coefficient_);
    struct FieldAlgoBaseInitData init_data(this->field_name_, 0, this->unit_si_, this->limits_, FieldFlag::Flags() );
    this->check_field_limits(this->in_rec_, init_data);

	return true;
}


// Instantion of Fields
INSTANCE_ALL(FieldTimeFunction)

