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
		.declare_key("time_function", TableFunction<Value>::get_input_type(), it::Default::obligatory(),
									"Values of time series initialization of Field.")
		.declare_key("unit", FieldAlgorithmBase<spacedim, Value>::get_input_type_unit_si(), it::Default::optional(),
									"Definition of unit.")
        .allow_auto_conversion("time_function")
		.close();
}

template <int spacedim, class Value>
const int FieldTimeFunction<spacedim, Value>::registrar =
		Input::register_class< FieldTimeFunction<spacedim, Value>, unsigned int >("FieldTimeFunction") +
		FieldTimeFunction<spacedim, Value>::get_input_type().size();


template <int spacedim, class Value>
FieldTimeFunction<spacedim, Value>::FieldTimeFunction( unsigned int n_comp)
: FieldConstant<spacedim, Value>(n_comp)
{}


template <int spacedim, class Value>
bool FieldTimeFunction<spacedim, Value>::set_time(const TimeStep &time)
{
	TableFunction<Value> table_function;

	if (!Value::is_scalable()) {
		WarningOut().fmt("Setting key 'time_function' of non-floating point field at address {}\nValues will be skipped.\n",
				this->in_rec_.address_string());
	}
	Input::Record func_rec = this->in_rec_;
	table_function.init_from_input( func_rec.val<Input::Record>("time_function") );
	this->r_value_ = table_function.value( time.end() );
    this->value_.scale(this->unit_conversion_coefficient_);

	return true;
}


// Instantion of Fields
INSTANCE_ALL(FieldTimeFunction)

