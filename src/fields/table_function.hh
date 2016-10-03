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
 * @file    table_function.hh
 * @brief
 */

#ifndef TABLE_FUNCTION_HH_
#define TABLE_FUNCTION_HH_

#include "input/input_type_forward.hh"
#include "input/input_exception.hh"
#include <vector>


/// Declaration of exceptions.
TYPEDEF_ERR_INFO( EI_LastTime, double);
TYPEDEF_ERR_INFO( EI_Time, double);
TYPEDEF_ERR_INFO( EI_MinTime, double);
TYPEDEF_ERR_INFO( EI_MaxTime, double);
DECLARE_INPUT_EXCEPTION( ExcNonAscendingTime, << "Nonascending order of declared time stamps in TableFunction:\n"
											  << EI_LastTime::val << " is followed by " << EI_Time::val << "." );
DECLARE_INPUT_EXCEPTION( ExcTimeOutOfRange,   << "Time " << EI_Time::val << " is out of range of TableFunction: <"
											  << EI_MinTime::val << ", " << EI_MaxTime::val << ">." );


template <class Value>
class TableFunction
{
public:
	/// Store value in one time stamp.
	struct TimeValue {
		TimeValue(double time, Value val)
		: time_(time), value_(r_value_), r_value_(val) {}

		double time_;
		Value value_;
		typename Value::return_type r_value_;
	};

    /**
     * Return record of one time stamp.
     */
    static const Input::Type::Tuple & get_input_type_time_val();

    /**
     * Return Record for time series initialization of Fields. Allow to use interpolation
     * of Field values defined in one field descriptor.
     */
    static const Input::Type::Record & get_input_type();

    /// Default constructor
    TableFunction();

    /// Initialize actual values of the field given from the given Input::Record @p rec.
    void init_from_input(const Input::Record &rec);

    /// Return true if TableFunction is initialized (method init_from_input was called).
    bool initialized();

    /**
     * Returns one value in one given point. ResultType can be used to avoid some costly calculation if the result is trivial.
     */
    typename Value::return_type const &value(double time);

private:
    /// Vector of values in all time stamps.
    std::vector<struct TimeValue> time_values_;

    /// Last time of computed value_ (to prevent repetitive calculation)
    double last_time_;

    /// Last value, prevents passing large values (vectors) by value.
    Value value_;
    typename Value::return_type r_value_;
};

#endif /* TABLE_FUNCTION_HH_ */
