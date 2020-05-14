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
 * @file    field_divide.hh
 * @brief
 */

#ifndef FIELD_DIVIDE_HH_
#define FIELD_DIVIDE_HH_

#include "fields/field.hh"


/**
 * This field is meant to be used to implement for initialization of velocity fields in Darcy flows.
 * You can either use velocity which is proportion of flux and cross section etc.
 * Unfortunately it introduce one more level of indirection, namely one more virtual call for getting the field value.
 *
 * - The field can not be initialized form the input.
 * - We allow only Scalar Value with element_type double of divisor field.
 */
template <int spacedim, class Value>
class FieldDivide : public FieldAlgorithmBase<spacedim, Value> {
public:
	typedef typename Field<spacedim, Value>::FactoryBase FactoryBaseType;
    typedef typename Space<spacedim>::Point Point;
    /**
     *
     */
    FieldDivide( std::shared_ptr< FieldAlgorithmBase<spacedim, Value> > inner_dividend, Field<3, FieldValue<3>::Scalar> inner_divisor, unsigned int n_comp=0);

    /**
     * Returns one value in one given point. ResultType can be used to avoid some costly calculation if the result is trivial.
     */
    virtual typename Value::return_type const &value(const Point &p, const ElementAccessor<spacedim> &elm);

    /**
     * Returns std::vector of scalar values in several points at once.
     */
    virtual void value_list (const Armor::array &point_list, const ElementAccessor<spacedim> &elm,
                       std::vector<typename Value::return_type>  &value_list);

    /**
     * Update time and possibly update data.
     */
    bool set_time(const TimeStep &time) override;

    virtual ~FieldDivide();

private:
    /// Field to which we add linear potential.
    std::shared_ptr< FieldAlgorithmBase<spacedim, Value> > inner_dividend_;
    /// Field to which we add linear potential.
    Field<3, FieldValue<3>::Scalar> inner_divisor_;
};


#endif /* FIELD_DIVIDE_HH_ */
