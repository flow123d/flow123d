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
 * @file    table_function.impl.hh
 * @brief
 */

#ifndef TABLE_FUNCTION_IMPL_HH_
#define TABLE_FUNCTION_IMPL_HH_

#include "fields/table_function.hh"
#include "input/input_type.hh"

namespace it = Input::Type;

template <class Value>
const it::Tuple & TableFunction<Value>::get_input_type_val()
{
    return it::Tuple("IndependentValue", "Value of Field for independent variable.")
        .declare_key("t", it::Double( 0.0 ), it::Default::obligatory(),
                                    "Independent variable of stamp." )
		.declare_key("value", Value::get_input_type(), it::Default::obligatory(),
									"Value of the field in given stamp." )
		.close();
}

template <class Value>
const it::Record & TableFunction<Value>::get_input_type()
{
    return it::Record("TableFunction", "Allow set variable series initialization of Fields.")
        .declare_key("values", it::Array( TableFunction<Value>::get_input_type_val(), 2 ), it::Default::obligatory(),
                                    "Initizaliation values of Field." )
        .allow_auto_conversion("values")
		.close();

}


template <class Value>
TableFunction<Value>::TableFunction()
: last_time_(-1.0),
  value_(r_value_)
{}


template <class Value>
void TableFunction<Value>::init_from_input(const Input::Record &rec)
{
	ASSERT( !this->initialized() ).error("TableFunction can't be initialized more than once.");

	Input::Array data_array = rec.val<Input::Array>("values");
	double last_time = -1.0;
    for (Input::Iterator<Input::Tuple> it = data_array.begin<Input::Tuple>(); it != data_array.end(); ++it) {
    	double time = it->val<double>("t");
    	if (last_time >= time) {
    		THROW( ExcNonAscendingTime() << EI_LastTime(last_time) << EI_Time(time) );
    	}
    	last_time = time;

    	typename Value::return_type r_value;
    	Value value(r_value);
    	value.init_from_input( it->val<typename Value::AccessType>("value") );
    	time_values_.push_back( TimeValue(time, value) );
    }
}

template <class Value>
bool TableFunction<Value>::initialized()
{
	return ( time_values_.size() > 0 );
}

template <class Value>
typename Value::return_type const &TableFunction<Value>::value(double time)
{
	ASSERT( this->initialized() ).error("Compute value of uninitialized TableFunction.");

	if (time != last_time_) {
		unsigned int last_idx = time_values_.size() - 1;
		if (time < time_values_[0].time_ || time > time_values_[last_idx].time_) {
			THROW( ExcTimeOutOfRange() << EI_Time(time) << EI_MinTime(time_values_[0].time_) << EI_MaxTime(time_values_[last_idx].time_) );
		}
		for (unsigned int i=0; i<last_idx; ++i) {
			if (time >= time_values_[i].time_ && time <= time_values_[i+1].time_) {
				double coef = (time - time_values_[i].time_) / (time_values_[i+1].time_ - time_values_[i].time_);
				this->value_.interpolated(coef, time_values_[i].r_value_, time_values_[i+1].r_value_);
				break;
			}
		}
	}

	return this->r_value_;
}

#endif /* TABLE_FUNCTION_IMPL_HH_ */
