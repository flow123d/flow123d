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

#include <string.h>                        // for memcpy
#include <vector>                          // for vector
#include <armadillo>
#include "fields/field_values.hh"          // for FieldValue<>::TensorFixed

namespace Input {
	class Record;
	namespace Type {
		class Record;
		class Tuple;
	}
}
class TimeStep;



template <class Value>
class TableFunction
{
public:
	typedef typename Value::return_type return_type;

	/// Store value in one t stamp.
	struct TableValue {
		TableValue(double t, Value val)
		: t_(t), value_(r_value_), r_value_(val) {}

		double t_;
		Value value_;
		return_type r_value_;
	};

    /**
     * Return record of one t stamp.
     */
    static const Input::Type::Tuple & get_input_type_val();

    /**
     * Return Record for one t stamp series initialization of Fields. Allow to use interpolation
     * of Field values defined in one field descriptor.
     */
    static const Input::Type::Record & get_input_type();

    /// Default constructor
    TableFunction();

    /// Initialize actual values of the field given from the given Input::Record @p rec.
    void init_from_input(const Input::Record &rec, const TimeStep &time);

    /// Return true if TableFunction is initialized (method init_from_input was called).
    bool initialized();

    /**
     * Returns one value in one given point. ResultType can be used to avoid some costly calculation if the result is trivial.
     */
    return_type const &value(double t);

private:
    // Compute the weighted average of table_values_[idx] and table_values_[idx+1]
    void interpolated(double coef, unsigned int idx);

    /// Vector of values in all stamps.
    std::vector<struct TableValue> table_values_;

    /// Last t stamp of computed value_ (to prevent repetitive calculation)
    double last_t_;

    /// Last value, prevents passing large values (vectors) by value.
    Value value_;
    return_type r_value_;
};

#endif /* TABLE_FUNCTION_HH_ */
