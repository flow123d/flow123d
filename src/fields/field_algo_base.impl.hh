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
#include "fields/field_interpolated_p0.hh"
#include "fields/field_python.hh"
#include "fields/field_constant.hh"
#include "fields/field_formula.hh"
#include "fields/field_elementwise.hh"

#include "fields/field_values.hh"

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
  component_idx_(std::numeric_limits<unsigned int>::max())
{
    value_.set_n_comp(n_comp);
}



template <int spacedim, class Value>
string FieldAlgorithmBase<spacedim, Value>::template_name() {
	return boost::str(boost::format("R%i -> %s") % spacedim % Value::type_name() );
}



template <int spacedim, class Value>
Input::Type::Abstract & FieldAlgorithmBase<spacedim, Value>::get_input_type() {
	stringstream ss;
	ss << "[" << Value::NRows_  << ", " << Value::NCols_  << "]";
    return it::Abstract("Field:"+template_name(), "Abstract for all time-space functions.")
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
shared_ptr< FieldAlgorithmBase<spacedim, Value> >
FieldAlgorithmBase<spacedim, Value>::function_factory(const Input::AbstractRecord &rec, unsigned int n_comp )
{
    shared_ptr< FieldAlgorithmBase<spacedim, Value> > func;
    func = rec.factory< FieldAlgorithmBase<spacedim, Value> >(n_comp);
    func->init_from_input(rec);
    return func;
}



template <int spacedim, class Value>
void FieldAlgorithmBase<spacedim, Value>::init_from_input(const Input::Record &rec) {
    xprintf(PrgErr, "The field '%s' do not support initialization from input.\n",
            typeid(this).name());
}



template <int spacedim, class Value>
bool FieldAlgorithmBase<spacedim, Value>::set_time(const TimeStep &time) {
    time_ = time;
    return false; // no change
}



template <int spacedim, class Value>
void FieldAlgorithmBase<spacedim, Value>::set_mesh(const Mesh *mesh,  bool boundary_domain) {
}



template<int spacedim, class Value>
unsigned int FieldAlgorithmBase<spacedim, Value>::n_comp() const {
    return (Value::NRows_ ? 0 : value_.n_rows());
}



template<int spacedim, class Value>
void FieldAlgorithmBase<spacedim, Value>::value_list(
        const std::vector< Point >  &point_list,
        const ElementAccessor<spacedim> &elm,
        std::vector<typename Value::return_type>  &value_list)
{
	ASSERT_EQ( point_list.size(), value_list.size() ).error();
    for(unsigned int i=0; i< point_list.size(); i++) {
    	ASSERT( Value(value_list[i]).n_rows()==this->value_.n_rows() )(i)(Value(value_list[i]).n_rows())(this->value_.n_rows())
                .error("value_list has wrong number of rows");
        value_list[i]=this->value(point_list[i], elm);
    }

}




#endif //FUNCTION_BASE_IMPL_HH_
