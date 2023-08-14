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
 * @file    table_function.cc
 * @brief
 */

#include "fields/table_function.hh"
#include "tools/time_governor.hh"
#include "input/input_type.hh"
#include "system/logger.hh"

namespace it = Input::Type;

template <class Value>
const it::Tuple & TableFunction<Value>::get_input_type_val()
{
	/*
	 * Table function is now fixed for representation of time value. It can be changed
	 * to independent value. It needs only replace type of 't' with generic type
	 * (Input::Type::Parameter).
	 */
    typedef FieldValue<3>::TensorFixed TensorValue;
    return it::Tuple("IndependentValue", "Value of Field for time variable.")
                                       //"Value of Field for independent variable."
        .declare_key("t", TimeGovernor::get_input_time_type( 0.0 ), it::Default::obligatory(),
                                    "Time stamp." )
                                  //"Independent variable of stamp."
		.declare_key("value", TensorValue::get_input_type(), it::Default::obligatory(),
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
: last_t_(-1.0),
  value_(r_value_)
{}


template <class Value>
void TableFunction<Value>::init_from_input(const Input::Record &rec, const TimeStep &time)
{
	ASSERT( !this->initialized() ).error("TableFunction can't be initialized more than once.");

	Input::Array data_array = rec.val<Input::Array>("values");
	double last_t = -1.0;
    for (Input::Iterator<Input::Tuple> it = data_array.begin<Input::Tuple>(); it != data_array.end(); ++it) {
    	double t = time.read_time(it->find<Input::Tuple>("t"));
    	if (last_t >= t) {
    		WarningOut().fmt("Nonascending order of declared stamps in TableFunction at address {}.\nStamp {} will be skipped.",
    				rec.address_string(), t);
    	} else {
        	last_t = t;

        	typename Value::return_type r_value;
        	Value value(r_value);
        	value.init_from_input( it->val<Input::Array>("value") );
        	table_values_.push_back( TableValue(t, value) );
    	}
    }
}

template <class Value>
bool TableFunction<Value>::initialized()
{
	return ( table_values_.size() > 0 );
}

template <class Value>
typename TableFunction<Value>::return_type const &TableFunction<Value>::value(double t)
{
	ASSERT( this->initialized() ).error("Compute value of uninitialized TableFunction.");

	if (t != last_t_) {
		unsigned int last_idx = table_values_.size() - 1;
		if (t < table_values_[0].t_) {
    		WarningOut().fmt("Value of stamp {} is out of range of TableFunction: <{}, {}>. Extrapolation of minimal value will be used.",
    				t, table_values_[0].t_, table_values_[last_idx].t_);
    		this->interpolated(0.0, 0);
		} else if (t > table_values_[last_idx].t_) {
    		WarningOut().fmt("Value of stamp {} is out of range of TableFunction: <{}, {}>. Extrapolation of maximal value will be used.",
    				t, table_values_[0].t_, table_values_[last_idx].t_);
    		this->interpolated(1.0, last_idx-1);
		} else {
			for (unsigned int i=0; i<last_idx; ++i) {
				if (t >= table_values_[i].t_ && t <= table_values_[i+1].t_) {
					double coef = (t - table_values_[i].t_) / (table_values_[i+1].t_ - table_values_[i].t_);
					this->interpolated(coef, i);
					break;
				}
			}
		}
	}

	return this->r_value_;
}

template <class Value>
void TableFunction<Value>::interpolated(double coef, unsigned int idx)
{
	ASSERT(coef >= 0 && coef <= 1)(coef).error();
	ASSERT(idx >= 0 && idx <= table_values_.size()-2)(idx).error();

	Value val_0(table_values_[idx].r_value_);
	Value val_1(table_values_[idx+1].r_value_);
    if (Value::is_scalable())
        for( unsigned int row=0; row<value_.n_rows(); row++)
            for( unsigned int col=0; col<value_.n_cols(); col++) {
            	value_(row,col) = val_0(row,col) + coef * (val_1(row,col) - val_0(row,col));
            }
}


// Instantiation of TableFunction class
template class TableFunction<FieldValue<0>::Enum >;
template class TableFunction<FieldValue<0>::Integer >;
template class TableFunction<FieldValue<0>::Scalar >;
template class TableFunction<FieldValue<0>::Vector >; // temporary solution for computing more fields at once in python
template class TableFunction<FieldValue<2>::VectorFixed >;
template class TableFunction<FieldValue<2>::TensorFixed >;
template class TableFunction<FieldValue<3>::VectorFixed >;
template class TableFunction<FieldValue<3>::TensorFixed >;
