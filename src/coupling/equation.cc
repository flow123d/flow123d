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
 * @file    equation.cc
 * @brief   Abstract base class for equation clasess.
 * @author  Jan Brezina
 */

#include <petscmat.h>
#include "tools/time_governor.hh"


#include "equation.hh"
#include "system/system.hh"
#include "input/accessors.hh"
#include "fields/field_set.hh"
#include "fields/field_common.hh"
#include "fields/bc_field.hh"
#include "tools/unit_converter.hh"
#include "tools/unit_si.hh"




/*****************************************************************************************
 * Implementation of EqBase
 */

Input::Type::Record & EquationBase::record_template() {
    return Input::Type::Record("EquationBase_AUX", "Auxiliary record with keys common for equations. Should not be used.")
        .declare_key("time", TimeGovernor::get_input_type(), Input::Type::Default("{}"),
                    "Time governor setting.")
		    .close();
}

Input::Type::Record & EquationBase::user_fields_template(std::string equation_name) {
    return Input::Type::Record("EquationBase_user_field_AUX", "Auxiliary record with common key user_field. Should not be used.")
        .declare_key("user_fields", Input::Type::Array(
                    FieldSet::make_user_field_type(equation_name)),
                    Input::Type::Default::optional(),
                    "Input fields of the equation defined by user.")
	    .close();
}

EquationBase::EquationBase()
: equation_empty_(true),
  mesh_(NULL),
  time_(NULL),
  input_record_(),
  eq_fieldset_(nullptr)
{}



EquationBase::EquationBase(Mesh &mesh, const  Input::Record in_rec)
: equation_empty_(false),
  mesh_(&mesh),
  time_(NULL),
  input_record_(in_rec),
  eq_fieldset_(nullptr)
{}


void EquationBase::set_time_governor(TimeGovernor &time)
{
  time_ = &time;
}

double EquationBase::solved_time()
{
    return time_->t();
}

void EquationBase::init_user_fields(Input::Array user_fields, FieldSet &output_fields) {
	for (Input::Iterator<Input::Record> it = user_fields.begin<Input::Record>();
                    it != user_fields.end();
                    ++it) {
	    std::string field_name = it->val<std::string>("name");
    	bool is_bdr = it->val<bool>("is_boundary");
    	auto shape_type = it->val<FieldSet::UserFieldShape>("shape_type");

    	// check if field of same name doesn't exist in FieldSet
    	auto * exist_field = eq_fieldset_->field(field_name);
    	if (exist_field!=nullptr) {
    	    THROW(FieldSet::ExcFieldExists() << FieldCommon::EI_Field(field_name));
    	}

    	UnitSI units = UnitSI::dimensionless();
    	Input::Record unit_record;
        if ( it->opt_val("unit", unit_record) ) {
            std::string unit_str = unit_record.val<std::string>("unit_formula");
        	try {
        		units.convert_unit_from(unit_str);
        	} catch (ExcInvalidUnit &e) {
        		e << it->ei_address();
        		throw;
        	} catch (ExcNoncorrespondingUnit &e) {
        		e << it->ei_address();
        		throw;
        	}
        }

    	Input::AbstractRecord field_rec = it->val<Input::AbstractRecord>("field");
        switch (shape_type)
        {
        case FieldSet::scalar:
            Field<3, FieldValue<3>::Scalar> * scalar_field;
            if (is_bdr)
                scalar_field = new BCField<3, FieldValue<3>::Scalar>();
            else
                scalar_field = new Field<3, FieldValue<3>::Scalar>();
            *eq_fieldset_+=scalar_field
                    ->name(field_name)
                    .description("")
                    .units( units )
					.flags(FieldFlag::equation_result);
            scalar_field->set_mesh(*mesh_);
            scalar_field->set( field_rec, time_->t());
            scalar_field->set_default_fieldset(*eq_fieldset_);
            output_fields+=*scalar_field;
            break;
        case FieldSet::vector:
            Field<3, FieldValue<3>::VectorFixed> * vector_field;
            if (is_bdr)
                vector_field = new BCField<3, FieldValue<3>::VectorFixed>();
            else
                vector_field = new Field<3, FieldValue<3>::VectorFixed>();
            *eq_fieldset_+=vector_field
                    ->name(field_name)
                    .description("")
                    .units( units )
					.flags(FieldFlag::equation_result);
            vector_field->set_mesh(*mesh_);
            vector_field->set( field_rec, time_->t());
            vector_field->set_default_fieldset(*eq_fieldset_);
            output_fields+=*vector_field;
            break;
        case FieldSet::tensor:
            Field<3, FieldValue<3>::TensorFixed> * tensor_field;
            if (is_bdr)
                tensor_field = new BCField<3, FieldValue<3>::TensorFixed>();
            else
                tensor_field = new Field<3, FieldValue<3>::TensorFixed>();
            *eq_fieldset_+=tensor_field
                    ->name(field_name)
                    .description("")
                    .units( units )
					.flags(FieldFlag::equation_result);
            tensor_field->set_mesh(*mesh_);
            tensor_field->set( field_rec, time_->t());
            tensor_field->set_default_fieldset(*eq_fieldset_);
            output_fields+=*tensor_field;
            break;
        }
	}
}
