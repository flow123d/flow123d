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
 * @file    multi_field.impl.hh
 * @brief   
 */

#ifndef MULTI_FIELD_IMPL_HH_
#define MULTI_FIELD_IMPL_HH_


#include "multi_field.hh"
#include "fields/field_algo_base.hh"
#include "fem/integral_acc.hh"
#include "fem/eval_points.hh"
#include "fem/element_cache_map.hh"
#include "fields/field_set.hh"
#include "input/input_exception.hh"
#include "io/observe.hh"

namespace it = Input::Type;

/******************************************************************************************
 * Implementation of MultiField<...>
 */

template<int spacedim, class Value>
MultiField<spacedim, Value>::MultiField()
: FieldCommon(),
  no_check_control_field_(nullptr)
{
// 	static_assert(Value::NRows_ == 1 && Value::NCols_ == 1, "");
	this->multifield_ = true;
    this->set_shape( Value::NRows_, Value::NCols_ );
}



template<int spacedim, class Value>
MultiField<spacedim, Value>::MultiField(const MultiField &other)
: FieldCommon(other),
  sub_fields_(other.sub_fields_),
  full_input_list_(other.full_input_list_),
  no_check_control_field_(other.no_check_control_field_)
{
	this->multifield_ = true;
	this->set_shape( Value::NRows_, Value::NCols_ );
}



template<int spacedim, class Value>
MultiField<spacedim,Value> &MultiField<spacedim,Value>::operator=(const MultiField<spacedim,Value> &other)
{
	//ASSERT( flags().match( FieldFlag::input_copy )(this->name())(other.name()).error("Try to assign to non-copy field.");
	ASSERT_PTR(other.shared_->mesh_).error("Must call set_mesh before assign to other field.\n");
	ASSERT( !shared_->mesh_ || (shared_->mesh_==other.shared_->mesh_))
	        .error("Assignment between multi fields with different meshes.\n");
	ASSERT_GT( shared_->comp_names_.size(), 0).error("Vector of component names can't be empty!\n");
	ASSERT_EQ( shared_->comp_names_.size(), other.shared_->comp_names_.size())
	        .error("Both multi fields must have same size of vectors of component names.\n");

	// check for self assignement
	if (&other == this) return *this;

	// class members derived from FieldCommon
	std::vector< std::string > comp_names = shared_->comp_names_; // keep component names
	shared_ = other.shared_;
	shared_->comp_names_ = comp_names;
    shared_->is_fully_initialized_ = false;
	set_time_result_ = other.set_time_result_;
	last_time_ = other.last_time_;
	last_limit_side_ = other.last_limit_side_;
	is_jump_time_ = other.is_jump_time_;
	component_index_ = other.component_index_;
	this->multifield_ = true;
	this->shape_ = other.shape_;

	// class members of Field class
	if ( size() == 0 ) {
		// move subfields from other, set correct names
		sub_fields_.clear();
		sub_fields_.reserve( other.size() );
		for (unsigned int i=0; i<other.size(); ++i) {
			sub_fields_.push_back( other.sub_fields_[i] );
			if (this->shared_->comp_names_[i].length() == 0)
				THROW( Input::ExcInputMessage() << EI_Message("The field " + this->input_name()
																	+ " has set empty name of component.") );
			else {
				sub_fields_[i].name_ = this->shared_->comp_names_[i] + "_" + name();
			}
		}
	} else {
		THROW( ExcMessage() << EI_Message("Internal error. Assignment operator can't be used after call setup_component() method.") );
	}
	full_input_list_ = other.full_input_list_;
	no_check_control_field_ = other.no_check_control_field_;

	return *this;
}



template<int spacedim, class Value>
it::Instance MultiField<spacedim,Value>::get_input_type() {
	ASSERT_PERMANENT(false).error("This method can't be used for MultiField");

	it::Abstract abstract = it::Abstract();
	it::Instance inst = it::Instance( abstract, std::vector<it::TypeBase::ParameterPair>() );
	return inst;
}


template<int spacedim, class Value>
it::Array MultiField<spacedim,Value>::get_multifield_input_type() {
	it::Array type = it::Array( SubFieldBaseType::get_input_type_instance(shared_->input_element_selection_), 1);
	return type;
}


template<int spacedim, class Value>
auto MultiField<spacedim, Value>::disable_where(
        const MultiField<spacedim, typename FieldValue<spacedim>::Enum > &control_field,
        const vector<FieldEnum> &value_list) -> MultiField &
{
    no_check_control_field_=&control_field;
    shared_->no_check_values_=value_list;
    return *this;
}


template<int spacedim, class Value>
bool MultiField<spacedim, Value>::set_time(
		const TimeStep &time, LimitSide limit_side)
{
	// initialization of Multifield for first call
	if (sub_fields_.size() == 0) {
	    setup_components();
	}

	// set time for sub fields
	set_time_result_ = TimeStatus::constant;
	is_jump_time_=false;
	for( SubFieldType &field : sub_fields_) {
            if (field.set_time(time, limit_side))
                set_time_result_ = TimeStatus::changed;
            is_jump_time_ = is_jump_time_ ||  field.is_jump_time();
	}
    return (set_time_result_ == TimeStatus::changed);
}



template<int spacedim, class Value>
void MultiField<spacedim, Value>::set_mesh(const Mesh &mesh) {
	// test if mesh is not set
	if (shared_->mesh_ && shared_->mesh_ != &mesh) {
		THROW(ExcFieldMeshDifference() << EI_Field(name()) );
	}

    shared_->mesh_ = &mesh;
}


template<int spacedim, class Value>
void MultiField<spacedim, Value>::copy_from(const FieldCommon & other) {
    ASSERT( flags().match(FieldFlag::equation_input))(other.name())(this->name())
            .error("Can not copy to the non-copy field.");

    // do not use copy if the field have its own input
    if ( flags().match(FieldFlag::declare_input)
         && this->shared_->input_list_.size() != 0 ) return;

    if (typeid(other) == typeid(*this)) {
		auto  const &other_field = dynamic_cast<  MultiField<spacedim, Value> const &>(other);
		this->operator=(other_field);
	} else if (typeid(other) == typeid(SubFieldType)) {
		auto  const &other_field = dynamic_cast<  SubFieldType const &>(other);
		sub_fields_.resize(1);
		sub_fields_[0] = other_field;
	}
}



template<int spacedim, class Value>
void MultiField<spacedim, Value>::field_output(std::shared_ptr<OutputTime> stream, OutputTime::DiscreteSpace type)
{
    ASSERT_LT( type, OutputTime::N_DISCRETE_SPACES ).error();

    for (unsigned long index=0; index < this->size(); index++) {
        sub_fields_[index].compute_field_data( type, stream );
    }
}




template<int spacedim, class Value>
bool MultiField<spacedim, Value>::is_constant(Region reg) {
	bool const_all=true;
	for(auto &field : sub_fields_) const_all = const_all && field.is_constant(reg);
	return const_all;
}

template<int spacedim, class Value>
FieldResult MultiField<spacedim, Value>::field_result( RegionSet region_set) const
{
    ASSERT(true).error("Not used yet. Test it.");

    FieldResult result_all = result_none;
    for(auto &field : sub_fields_) {
        FieldResult sub_result = field.field_result(region_set);
        if (sub_result == result_none) return result_none;

        if (result_all == result_none) // first subfield
            result_all = sub_result;
        else if (sub_result == result_other || result_all == result_other)
            result_all = result_other;
        else if (sub_result != result_all)
            result_all = result_constant; // all subfields are (possibly different) constants
    }

    return result_all;

}


template<int spacedim, class Value>
std::string MultiField<spacedim, Value>::get_value_attribute() const
{
    int nrows = Value::NRows_;
    int ncols = Value::NCols_;
    string type = "Integer";
    if (std::is_floating_point<typename Value::element_type>::value)
        type = "Double";

    return fmt::format("{{ \"subfields\": true, \"shape\": [ {}, {} ], \"type\": \"{}\", \"limit\": [ {}, {} ] }}",
    					nrows, ncols, type, this->limits().first, this->limits().second);
}


template<int spacedim, class Value>
void MultiField<spacedim, Value>::setup_components() {
	unsigned int comp_size = this->shared_->comp_names_.size();
	ASSERT_GT(comp_size, 0).error("Vector of component names is empty!\n");
	ASSERT_PTR(this->shared_->mesh_).error("Mesh is not set!\n");

    sub_fields_.reserve( comp_size );
    for(unsigned int i_comp=0; i_comp < comp_size; i_comp++)
    {
    	sub_fields_.push_back( SubFieldType(i_comp, name(), this->full_comp_name(i_comp)) );
    	sub_fields_[i_comp].units( units() );
        if (no_check_control_field_ != nullptr && no_check_control_field_->size() == sub_fields_.size())
          sub_fields_[i_comp].disable_where((*no_check_control_field_)[i_comp], shared_->no_check_values_);
    	sub_fields_[i_comp].set_mesh( *(shared_->mesh_) );
//     	sub_fields_[i_comp].set_limit_side(this->limit_side_);
        sub_fields_[i_comp].input_selection(shared_->input_element_selection_);
    	sub_fields_[i_comp].add_factory( std::make_shared<MultiFieldFactory>(i_comp) );

    	if (this->shared_->input_default_!="") {
    		sub_fields_[i_comp].shared_->input_default_ = this->shared_->input_default_;
    	}

    	sub_fields_[i_comp].flags_ = this->flags_;
    	sub_fields_[i_comp].set_input_list(this->full_input_list_, *tg_);
    	sub_fields_[i_comp].set_default_fieldset( *(this->shared_->default_fieldset_) );
    }
}



template<int spacedim, class Value>
void MultiField<spacedim,Value>::set_input_list(const Input::Array &list, const TimeGovernor &tg) {
    if (! flags().match(FieldFlag::declare_input)) return;

    // Check sizes of Arrays defined MultiField in field descriptors
    for (Input::Iterator<Input::Record> it = list.begin<Input::Record>();
					it != list.end();
					++it) {
    	Input::Array mf_array;
    	if ( it->opt_val(this->input_name(), mf_array) ) {
    		unsigned int comp_size = this->shared_->comp_names_.size();
    		if (mf_array.size() != 1 && mf_array.size() != comp_size)
    			THROW( Exc_InvalidMultiFieldSize() << EI_MultiFieldName(this->input_name())
    					<< EI_Size(mf_array.size()) << EI_ExpectedSize(comp_size) << list.ei_address() );
    	}
    }
    
    this->full_input_list_ = list;
    this->tg_ = &tg;
    
    // Save the full array for future use in FieldCommon::mark_input_times().
    list.copy_to(shared_->input_list_);
}



template<int spacedim, class Value>
typename Field<spacedim,Value>::FieldBasePtr MultiField<spacedim, Value>::MultiFieldFactory::create_field(Input::Record descriptor_rec, const FieldCommon &field) {
	Input::Array multifield_arr;
	if (descriptor_rec.opt_val(field.input_name(), multifield_arr))
	{
		//ASSERT_PERMANENT(multifield_arr.size() == 1 || multifield_arr.size() == field.n_comp() )(field.input_name())
		//		.error("Invalid size of Array defined for MultiField!\n";
		unsigned int position = 0;
		auto it = multifield_arr.begin<Input::AbstractRecord>();
		if (multifield_arr.size() > 1)
			while (index_ != position) {
				++it; ++position;
			}

		FieldAlgoBaseInitData init_data(field.input_name(), field.n_comp(), field.units(), field.limits(), field.get_flags());
		typename Field<spacedim,Value>::FieldBasePtr field_algo_base = Field<spacedim,Value>::FieldBaseType::function_factory( (*it), init_data );
		field_algo_base->set_component_idx(index_);
		return field_algo_base;
	}

	return NULL;
}



template<int spacedim, class Value>
bool MultiField<spacedim, Value>::MultiFieldFactory::is_active_field_descriptor(const Input::Record &in_rec, const std::string &input_name) {
	return in_rec.find<Input::Array>(input_name);
}


template<int spacedim, class Value>
std::vector<const FieldCommon *> MultiField<spacedim, Value>::set_dependency(FMT_UNUSED unsigned int i_reg) const {
    ASSERT_PERMANENT(false).error("Set dependency of MultiField should be performed by individual components!\n");
    return std::vector<const FieldCommon *>();
}



template<int spacedim, class Value>
void MultiField<spacedim, Value>::cache_reallocate(FMT_UNUSED const ElementCacheMap &cache_map, FMT_UNUSED unsigned int region_idx) const {
    //for(auto &field : sub_fields_) field.cache_reallocate(cache_map);
	ASSERT_PERMANENT(false).error("Cache reallocate of MultiField should be performed by individual components!\n");
}


template<int spacedim, class Value>
void MultiField<spacedim, Value>::cache_update(FMT_UNUSED ElementCacheMap &cache_map,
        FMT_UNUSED unsigned int region_patch_idx) const {
    ASSERT_PERMANENT(false).error("Cache update of MultiField should be performed by individual components!\n");
}


template<int spacedim, class Value>
void MultiField<spacedim, Value>::set(
        std::vector<typename Field<spacedim, Value>::FieldBasePtr> field_vec,
        double time,
		std::vector<std::string> region_set_names)
{
	unsigned int comp_size = this->shared_->comp_names_.size();
	ASSERT_GT(comp_size, 0).error("Vector of component names is empty!\n");
	ASSERT_EQ(comp_size, field_vec.size());
	ASSERT_PTR(this->shared_->mesh_).error("Mesh is not set!\n");

    sub_fields_.reserve( comp_size );
    for(unsigned int i_comp=0; i_comp < comp_size; i_comp++)
    {
    	sub_fields_.push_back( SubFieldType(i_comp, name(), "") );
    	sub_fields_[i_comp].set_mesh( *(shared_->mesh_) );
    	sub_fields_[i_comp].flags_ = this->flags_;
    	sub_fields_[i_comp].set(field_vec[i_comp], time, region_set_names);
    }
}


template<int spacedim, class Value>
void MultiField<spacedim, Value>::set(
        typename Field<spacedim, Value>::FieldBasePtr field,
        double time,
		std::vector<std::string> region_set_names)
{
	ASSERT_EQ(this->shared_->comp_names_.size(), 1).error("Size of component names vector must be 1!\n");
	ASSERT_PTR(this->shared_->mesh_).error("Mesh is not set!\n");

    sub_fields_.reserve(1);
   	sub_fields_.push_back( SubFieldType(0, name(), "") );
   	sub_fields_[0].set_mesh( *(shared_->mesh_) );
   	sub_fields_[0].flags_ = this->flags_;
   	sub_fields_[0].set(field, time, region_set_names);
}



#endif /* MULTI_FIELD_IMPL_HH_ */
