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
 * @file    field_divide.cc
 * @brief
 */

#include "fields/field_divide.hh"
#include "fields/field_instances.hh"	// for instantiation macros


template <int spacedim, class Value>
FieldDivide<spacedim, Value>::FieldDivide( std::shared_ptr< FieldAlgorithmBase<spacedim, Value> > inner_dividend,
        Field<3, FieldValue<3>::Scalar> inner_divisor, unsigned int n_comp)
: FieldAlgorithmBase<spacedim, Value>(n_comp),
  inner_dividend_(inner_dividend),
  inner_divisor_(inner_divisor)
{
    inner_divisor_.flags_add(FieldFlag::declare_input);
}



template <int spacedim, class Value>
typename Value::return_type const & FieldDivide<spacedim, Value>::value(const Point &p, const ElementAccessor<spacedim> &elm)
{
    this->r_value_ = inner_dividend_->value(p,elm);
    double div_val = inner_divisor_.value(p,elm);

    for(unsigned int row=0; row < this->value_.n_rows(); row++)
        for(unsigned int col=0; col < this->value_.n_cols(); col++) {
            this->value_(row,col) /= div_val;
        }

    return this->r_value_;
}



/**
 * Returns std::vector of scalar values in several points at once.
 */
template <int spacedim, class Value>
void FieldDivide<spacedim, Value>::value_list (const Armor::array &point_list, const ElementAccessor<spacedim> &elm,
                   std::vector<typename Value::return_type>  &value_list)
{
	ASSERT_EQ( point_list.size(), value_list.size() ).error();
    ASSERT_DBG( point_list.n_rows() == spacedim && point_list.n_cols() == 1 ).error("Invalid point size.\n");

	inner_dividend_->value_list(point_list, elm, value_list);

    for(unsigned int i=0; i< point_list.size(); i++) {
        double div_val = inner_divisor_.value(point_list.vec<spacedim>(i),elm);
        Value envelope(value_list[i]);

        for(unsigned int row=0; row < this->value_.n_rows(); row++)
            for(unsigned int col=0; col < this->value_.n_cols(); col++)
                envelope(row,col) /= div_val;
    }
}



template <int spacedim, class Value>
bool FieldDivide<spacedim,Value>::set_time (const TimeStep &time)
{
	ASSERT_PTR(inner_dividend_).error("Null data pointer.\n");
	bool changed = inner_divisor_.set_time(time, LimitSide::right);
    return inner_dividend_->set_time(time) && changed;
}



template <int spacedim, class Value>
FieldDivide<spacedim, Value>::~FieldDivide() {
}


// Instantiations of FieldConstant
INSTANCE_ALL(FieldDivide)

