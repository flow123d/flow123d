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
 * @file    field_algo_base.impl.hh
 * @brief   
 */

#ifndef field_algo_base_IMPL_HH_
#define field_algo_base_IMPL_HH_

#include <string>
#include <limits>
#include <memory>
using namespace std;

#include "fields/field_algo_base.hh"
#include "fields/field_python.hh"
#include "fields/field_constant.hh"
#include "fields/field_formula.hh"

#include "fields/field_values.hh"

#include "tools/unit_converter.hh"

#include "tools/time_governor.hh"
#include "input/factory.hh"
#include "input/accessors.hh"
#include "input/flow_attribute_lib.hh"

namespace it = Input::Type;




/******************************************************************************************
 * Implementation of FieldBase<...>
 */

template <int spacedim, class Value>
FieldAlgorithmBase<spacedim, Value>::FieldAlgorithmBase(unsigned int n_comp)
: value_(r_value_),
  field_result_(result_other),
  component_idx_(std::numeric_limits<unsigned int>::max()),
  unit_conversion_coefficient_(1.0)
{
    value_.set_n_comp(n_comp);
}



template <int spacedim, class Value>
string FieldAlgorithmBase<spacedim, Value>::template_name() {
	return boost::str(boost::format("R%i_to_%s") % spacedim % Value::type_name() );
}



template <int spacedim, class Value>
Input::Type::Abstract & FieldAlgorithmBase<spacedim, Value>::get_input_type() {
	stringstream ss;
	ss << "[" << Value::NRows_  << ", " << Value::NCols_  << "]";
    return it::Abstract("Field_"+template_name(), "Abstract for all time-space functions.")
			.allow_auto_conversion("FieldConstant")
			.root_of_generic_subtree()
			.add_attribute(FlowAttribute::field_value_shape(), ss.str() )
			.close();
}


template <int spacedim, class Value>
const Input::Type::Instance & FieldAlgorithmBase<spacedim, Value>::get_input_type_instance(Input::Type::Selection value_selection) {
	std::vector<it::TypeBase::ParameterPair> param_vec;
	if (is_enum_valued) {
		ASSERT( !(value_selection==Input::Type::Selection()) ).error("Not defined 'value_selection' for enum element type.\n");
		param_vec.push_back( std::make_pair("element_input_type", std::make_shared<it::Selection>(value_selection)) );
	} else {
		param_vec.push_back( std::make_pair("element_input_type", std::make_shared<typename Value::ElementInputType>()) );
	}

	return it::Instance(get_input_type(), param_vec).close();
}

template <int spacedim, class Value>
const Input::Type::Record & FieldAlgorithmBase<spacedim, Value>::get_field_algo_common_keys() {
    auto unit_record = it::Record("Unit",
           "Specify unit of an input value. "
           "Evaluation of the unit formula results into a coeficient and a "
           "unit in terms of powers of base SI units. The unit must match "
           "expected SI unit of the value, while the value provided on the input "
           "is multiplied by the coefficient before further processing. "
           "The unit formula have a form:\n"
           "```\n"
           "<UnitExpr>;<Variable>=<Number>*<UnitExpr>;...,\n"
           "```\n"
           "where ```<Variable>``` is a variable name and ```<UnitExpr>``` is a units expression "
           "which consists of products and divisions of terms.\n\n"
           "A term has a form: "
           "```<Base>^<N>```, where ```<N>``` is an integer exponent and ```<Base>``` "
           "is either a base SI unit, a derived unit, or a variable defined in the same unit formula. "
           "Example, unit for the pressure head:\n\n"
           "```MPa/rho/g_; rho = 990*kg*m^-3; g_ = 9.8*m*s^-2```"
            )
        .allow_auto_conversion("unit_formula")
        .declare_key("unit_formula", it::String(), it::Default::obligatory(),
                                   "Definition of unit." )
        .close();

    return it::Record("FieldAlgorithmBase_common_aux", "")
        .declare_key("unit", unit_record, it::Default::optional(),
                                "Unit of the field values provided in the main input file, in the external file, or "
                                "by a function (FieldPython).")
        .close();
}

template <int spacedim, class Value>
shared_ptr< FieldAlgorithmBase<spacedim, Value> >
FieldAlgorithmBase<spacedim, Value>::function_factory(const Input::AbstractRecord &rec, const struct FieldAlgoBaseInitData& init_data )
{
    shared_ptr< FieldAlgorithmBase<spacedim, Value> > func;
    func = rec.factory< FieldAlgorithmBase<spacedim, Value> >(init_data.n_comp_);
    func->init_from_input(rec, init_data);
    return func;
}



template <int spacedim, class Value>
void FieldAlgorithmBase<spacedim, Value>::init_from_input(const Input::Record &, const struct FieldAlgoBaseInitData&) {
    xprintf(PrgErr, "The field '%s' do not support initialization from input.\n",
            typeid(this).name());
}



template <int spacedim, class Value>
bool FieldAlgorithmBase<spacedim, Value>::set_time(const TimeStep &time) {
    time_ = time;
    return false; // no change
}



template <int spacedim, class Value>
void FieldAlgorithmBase<spacedim, Value>::set_mesh(const Mesh *,  bool) {
}



template<int spacedim, class Value>
unsigned int FieldAlgorithmBase<spacedim, Value>::n_comp() const {
    return (Value::NRows_ ? 0 : value_.n_rows());
}



template<int spacedim, class Value>
void FieldAlgorithmBase<spacedim, Value>::value_list(
        const Armor::array  &point_list,
        const ElementAccessor<spacedim> &elm,
        std::vector<typename Value::return_type>  &value_list)
{
	ASSERT_EQ( point_list.size(), value_list.size() ).error();
    ASSERT_DBG(point_list.n_rows() == spacedim && point_list.n_cols() == 1).error("Invalid point size.\n");
    for(unsigned int i=0; i< point_list.size(); i++) {
    	ASSERT( Value(value_list[i]).n_rows()==this->value_.n_rows() )(i)(Value(value_list[i]).n_rows())(this->value_.n_rows())
                .error("value_list has wrong number of rows");
        value_list[i]=this->value(point_list.vec<spacedim>(i), elm);
    }

}

template<int spacedim, class Value>
void FieldAlgorithmBase<spacedim, Value>::init_unit_conversion_coefficient(const Input::Record &rec,
		const struct FieldAlgoBaseInitData& init_data)
{
    Input::Record unit_record;
    if ( rec.opt_val("unit", unit_record) ) {
        if (!Value::is_scalable()) {
            WarningOut().fmt("Setting unit conversion coefficient of non-floating point field at address {}\nCoefficient will be skipped.\n",
                    rec.address_string());
        }
        std::string unit_str = unit_record.val<std::string>("unit_formula");
    	try {
    		this->unit_conversion_coefficient_ = init_data.unit_si_.convert_unit_from(unit_str);
    	} catch (ExcInvalidUnit &e) {
    		e << rec.ei_address();
    		throw;
    	} catch (ExcNoncorrespondingUnit &e) {
    		e << rec.ei_address();
    		throw;
    	}
    }
}




#endif //FUNCTION_BASE_IMPL_HH_
