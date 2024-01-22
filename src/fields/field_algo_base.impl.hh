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
  component_idx_(undef_idx),
  unit_conversion_coefficient_(1.0)
{
    value_.set_n_comp(n_comp);
}



template <int spacedim, class Value>
string FieldAlgorithmBase<spacedim, Value>::template_name() {
	return ""; //fmt::format("R{:d}_to_{}", spacedim, Value::type_name() ); // maybe change template name
}



template <int spacedim, class Value>
Input::Type::Abstract & FieldAlgorithmBase<spacedim, Value>::get_input_type() {
	//stringstream ss;
	//ss << "[" << Value::NRows_  << ", " << Value::NCols_  << "]";
    return it::Abstract("Field_"+template_name(), "Abstract for all time-space functions.")
			.allow_auto_conversion("FieldConstant")
			.root_of_generic_subtree()
			//.add_attribute(FlowAttribute::field_value_shape(), ss.str() )
			.close();
}


template <int spacedim, class Value>
const Input::Type::Instance & FieldAlgorithmBase<spacedim, Value>::get_input_type_instance(Input::Type::Selection value_selection) {
	std::vector<it::TypeBase::ParameterPair> param_vec;
	if (is_enum_valued) {
		ASSERT_PERMANENT( !(value_selection==Input::Type::Selection()) ).error("Not defined 'value_selection' for enum element type.\n");
		param_vec.push_back( std::make_pair("element_input_type", std::make_shared<it::Selection>(value_selection)) );
	} else {
		param_vec.push_back( std::make_pair("element_input_type", std::make_shared<typename Value::ElementInputType>()) );
	}

	return it::Instance(get_input_type(), param_vec).close();
}

template <int spacedim, class Value>
const Input::Type::Record & FieldAlgorithmBase<spacedim, Value>::get_field_algo_common_keys() {
    return it::Record("FieldAlgorithmBase_common_aux", "")
        .declare_key("unit", UnitConverter::get_input_type(), it::Default::optional(),
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
    THROW( ExcInputInitUnsupported() << EI_Field(typeid(this).name()) );
}



template <int spacedim, class Value>
bool FieldAlgorithmBase<spacedim, Value>::set_time(const TimeStep &time) {
    time_ = time;
    return false; // no change
}



template <int spacedim, class Value>
void FieldAlgorithmBase<spacedim, Value>::set_mesh(const Mesh *) {
}



template<int spacedim, class Value>
unsigned int FieldAlgorithmBase<spacedim, Value>::n_comp() const {
    return (Value::NRows_ ? 0 : value_.n_rows());
}


template<int spacedim, class Value>
void FieldAlgorithmBase<spacedim, Value>::cache_update(
            FMT_UNUSED FieldValueCache<typename Value::element_type> &data_cache,
			FMT_UNUSED ElementCacheMap &cache_map,
			FMT_UNUSED unsigned int region_patch_idx)
{
    ASSERT_PERMANENT(false).error("Must be implemented in descendants!\n");
}


template<int spacedim, class Value>
void FieldAlgorithmBase<spacedim, Value>::cache_reinit(FMT_UNUSED const ElementCacheMap &cache_map)
{}


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
