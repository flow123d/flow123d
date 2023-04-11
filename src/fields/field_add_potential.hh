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
 * @file    field_add_potential.hh
 * @brief   
 */

#ifndef FIELD_ADD_POTENTIAL_HH_
#define FIELD_ADD_POTENTIAL_HH_

#include <armadillo>
#include <memory>

#include "fields/field.hh"
#include "fields/field_algo_base.hh"
#include "fields/field_model.hh"
#include "fields/field_coords.hh"


/*******************************************************************************
 * Functors of FieldModels
 */
using Sclr = double;
using Vect = arma::vec3;

// Functor computing piezo_head_p0
struct fn_add_potential {
	inline Sclr operator() (Vect gravity, Vect coords, Sclr pressure) {
        return arma::dot(gravity, coords) + pressure;
    }
};


/**
 * Factory class (descendant of @p Field<...>::FactoryBase) that is necessary
 * for setting pressure values are piezometric head values.
 */
template <int spacedim, class Value>  // <3, FieldValue<3>::Scalar>
class AddPotentialFactory : public Field<spacedim, Value>::FactoryBase {
public:
    /// Constructor.
    AddPotentialFactory( Field<3, FieldValue<3>::VectorFixed > &gravity, FieldCoords &coords, Field<3, FieldValue<3>::Scalar> &inner_field)
    : gravity_(gravity),
	  coords_(coords),
	  inner_field_(inner_field),
      field_name_(inner_field.input_name())
    {}

    typename Field<spacedim,Value>::FieldBasePtr create_field(Input::Record rec, const FieldCommon &) override {
        Input::AbstractRecord field_a_rec;
        if (rec.opt_val(field_name_, field_a_rec)) {
           	return Model<3, FieldValue<3>::Scalar>::create(fn_add_potential(), gravity_, coords_, inner_field_);

        } else {
            return typename Field<spacedim,Value>::FieldBasePtr();
        }
    }

    bool is_active_field_descriptor(const Input::Record &in_rec, FMT_UNUSED const std::string &input_name) override {
        return in_rec.find<Input::AbstractRecord>(field_name_);
    }

    Field<3, FieldValue<3>::VectorFixed > &gravity_;
    FieldCoords &coords_;
    Field<3, FieldValue<3>::Scalar> &inner_field_;
    std::string field_name_;
};

#endif /* FIELD_ADD_POTENTIAL_HH_ */
