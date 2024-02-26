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
 * @file    field.impl.hh
 * @brief   
 */

#ifndef FIELD_IMPL_HH_
#define FIELD_IMPL_HH_

#include "field.hh"
#include "field_algo_base.impl.hh"
#include "field_fe.hh"
#include "fields/eval_subset.hh"
#include "fields/eval_points.hh"
#include "fields/field_value_cache.hh"
#include "fields/field_set.hh"
#include "mesh/region.hh"
#include "input/reader_to_storage.hh"
#include "input/accessors.hh"
#include "io/observe.hh"
#include "io/output_mesh.hh"
#include "io/output_element.hh"


/******************************************************************************************
 * Implementation of Field<...>
 */

template<int spacedim, class Value>
Field<spacedim,Value>::Field()
: data_(std::make_shared<SharedData>()),
  value_cache_( FieldValueCache<typename Value::element_type>(Value::NRows_, Value::NCols_) )
{
	// n_comp is nonzero only for variable size vectors Vector, VectorEnum, ..
	// this invariant is kept also by n_comp setter
	shared_->n_comp_ = (Value::NRows_ ? 0 : 1);
	this->add_factory( std::make_shared<FactoryBase>() );

	unsigned int cache_size = CacheMapElementNumber::get();
	value_cache_.reinit(cache_size);
	value_cache_.resize(cache_size);

	this->multifield_ = false;
	this->set_shape( Value::NRows_, Value::NCols_ );
}


template<int spacedim, class Value>
Field<spacedim,Value>::Field(const string &name)
: data_(std::make_shared<SharedData>()),
  value_cache_( FieldValueCache<typename Value::element_type>(Value::NRows_, Value::NCols_) )
{
		// n_comp is nonzero only for variable size vectors Vector, VectorEnum, ..
		// this invariant is kept also by n_comp setter
		shared_->n_comp_ = (Value::NRows_ ? 0 : 1);
		this->name( name );
		this->add_factory( std::make_shared<FactoryBase>() );
		unsigned int cache_size = CacheMapElementNumber::get();
		value_cache_.reinit(cache_size);
		value_cache_.resize(cache_size);

		this->multifield_ = false;
		this->set_shape( Value::NRows_, Value::NCols_ );
}



template<int spacedim, class Value>
Field<spacedim,Value>::Field(unsigned int component_index, string input_name, string name)
: data_(std::make_shared<SharedData>()),
  value_cache_( FieldValueCache<typename Value::element_type>(Value::NRows_, Value::NCols_) )
{
	// n_comp is nonzero only for variable size vectors Vector, VectorEnum, ..
	// this invariant is kept also by n_comp setter
	shared_->n_comp_ = (Value::NRows_ ? 0 : 1);
	this->set_component_index(component_index);
	this->name_ = (name=="") ? input_name : name;
	this->shared_->input_name_ = input_name;

	unsigned int cache_size = CacheMapElementNumber::get();
	value_cache_.reinit(cache_size);
	value_cache_.resize(cache_size);

	this->multifield_ = false;
	this->set_shape( Value::NRows_, Value::NCols_ );
}


template<int spacedim, class Value>
Field<spacedim,Value>::Field(const Field &other)
: FieldCommon(other),
  data_(other.data_),
  region_fields_(other.region_fields_),
  factories_(other.factories_),
  value_cache_(other.value_cache_)
{
	if (other.no_check_control_field_)
		no_check_control_field_ =  make_shared<ControlField>(*other.no_check_control_field_);

	this->multifield_ = false;
	this->set_shape( Value::NRows_, Value::NCols_ );
}


template<int spacedim, class Value>
Field<spacedim,Value> &Field<spacedim,Value>::operator=(const Field<spacedim,Value> &other)
{
	//ASSERT( flags().match(FieldFlag::input_copy) )(this->name())(other.name()).error("Try to assign to non-copy field from the field.");
	ASSERT(other.shared_->mesh_).error("Must call set_mesh before assign to other field.\n");
	ASSERT( !shared_->mesh_ || (shared_->mesh_==other.shared_->mesh_) ).error("Assignment between fields with different meshes.\n");

	// check for self assignement
	if (&other == this) return *this;

	// class members derived from FieldCommon
	shared_ = other.shared_;
    shared_->is_fully_initialized_ = false;
	set_time_result_ = other.set_time_result_;
	last_time_ = other.last_time_;
	last_limit_side_ = other.last_limit_side_;
	is_jump_time_ = other.is_jump_time_;
	component_index_ = other.component_index_;
	this->multifield_ = false;

	// class members of Field class
	data_ = other.data_;
	factories_ = other.factories_;
	region_fields_ = other.region_fields_;
	value_cache_ = other.value_cache_;
	this->shape_ = other.shape_;

	if (other.no_check_control_field_) {
		no_check_control_field_ =  make_shared<ControlField>(*other.no_check_control_field_);
	}

	return *this;
}



template<int spacedim, class Value>
typename Value::return_type Field<spacedim,Value>::operator() (BulkPoint &p) {
    return p.elm_cache_map()->get_value<Value>(value_cache_, p.elem_patch_idx(), p.eval_point_idx());
}



template<int spacedim, class Value>
typename Value::return_type Field<spacedim,Value>::operator() (SidePoint &p) {
    return p.elm_cache_map()->get_value<Value>(value_cache_, p.elem_patch_idx(), p.eval_point_idx());
}



template<int spacedim, class Value>
typename Value::return_type
Field<spacedim,Value>::operator[] (unsigned int i_cache_point) const
{
	return Value::get_from_array( this->value_cache_, i_cache_point );
}



template<int spacedim, class Value>
it::Instance Field<spacedim,Value>::get_input_type() {
	return FieldBaseType::get_input_type_instance(shared_->input_element_selection_);
}



template<int spacedim, class Value>
it::Array Field<spacedim,Value>::get_multifield_input_type() {
	ASSERT_PERMANENT(false).error("This method can't be used for Field");

	it::Array arr = it::Array( it::Integer() );
	return arr;
}



template<int spacedim, class Value>
auto Field<spacedim, Value>::disable_where(
		const Field<spacedim, typename FieldValue<spacedim>::Enum > &control_field,
		const vector<FieldEnum> &value_list) -> Field &
{
    no_check_control_field_=std::make_shared<ControlField>(control_field);
    shared_->no_check_values_=value_list;
    return *this;
}

template<int spacedim, class Value>
void Field<spacedim,Value>::set_mesh(const Mesh &in_mesh) {
	// since we allow copy of fields before set_mesh is called
	// we have to check that all copies set the same mesh
	if (shared_->mesh_ && shared_->mesh_ != &in_mesh) {
		THROW(ExcFieldMeshDifference() << EI_Field(name()) );
	}

	shared_->mesh_ = &in_mesh;

    // initialize table if it is empty, we assume that the RegionDB is closed at this moment
	region_fields_.resize( mesh()->region_db().size() );
	RegionHistory init_history(history_length_limit_);	// capacity
    data_->region_history_.resize( mesh()->region_db().size(), init_history );

    if (no_check_control_field_) no_check_control_field_->set_mesh(in_mesh);
}


/*
template<int spacedim, class Value>
std::shared_ptr< typename Field<spacedim,Value>::FieldBaseType >
Field<spacedim,Value>::operator[] (Region reg)
{
    ASSERT_PERMANENT_LT(reg.idx(), this->region_fields_.size());
    return this->region_fields_[reg.idx()];
}
*/


template <int spacedim, class Value>
bool Field<spacedim, Value>::is_constant(Region reg) {
	ASSERT(this->set_time_result_ != TimeStatus::unknown).error("Unknown time status.\n");
	ASSERT_LT(reg.idx(), this->region_fields_.size());
    FieldBasePtr region_field = this->region_fields_[reg.idx()];
    return ( region_field && region_field->is_constant_in_space() );
}


template<int spacedim, class Value>
void Field<spacedim, Value>::set(
		FieldBasePtr field,
		double time,
		std::vector<std::string> region_set_names)
{
	ASSERT_PTR(field).error("Null field pointer.\n");

	ASSERT_PTR( mesh() ).error("Null mesh pointer, set_mesh() has to be called before set().\n");
    ASSERT_EQ( field->n_comp() , shared_->n_comp_);
    field->set_mesh( mesh() );

	for (auto set_name : region_set_names) {
		RegionSet domain = mesh()->region_db().get_region_set(set_name);
        if (domain.size() == 0) continue;

        HistoryPoint hp = HistoryPoint(time, field);
        for(const Region &reg: domain) {
    	    RegionHistory &region_history = data_->region_history_[reg.idx()];
    	    // insert hp into descending time sequence
    	    ASSERT( region_history.size() == 0 || region_history[0].first < hp.first)(hp.first)(region_history[0].first)
    		        .error("Can not insert smaller time then last time in field's history.\n");
    	    region_history.push_front(hp);
        }
	}
    set_history_changed();
}



template<int spacedim, class Value>
void Field<spacedim, Value>::set(
        const Input::AbstractRecord &a_rec,
        double time,
		std::vector<std::string> region_set_names)
{
	FieldAlgoBaseInitData init_data(input_name(), n_comp(), units(), limits(), flags());
	set(FieldBaseType::function_factory(a_rec, init_data), time, region_set_names);
}



template<int spacedim, class Value>
bool Field<spacedim, Value>::set_time(const TimeStep &time_step, LimitSide limit_side)
{
	ASSERT_PTR( mesh() )( name() ).error("NULL mesh pointer of field with given name. set_mesh must be called before.\n");

    // Skip setting time if the new time is equal to current time of the field
	// and if either the field is continuous in that time or the current limit side is same as the new one.
    if (time_step.end() == last_time_) {
        if ( ! is_jump_time() ||
             limit_side == last_limit_side_) {
            last_limit_side_ = limit_side;
            return changed();
        }
    }

    last_time_=time_step.end();
    last_limit_side_ = limit_side;

    // possibly update our control field
    if (no_check_control_field_) {
            no_check_control_field_->set_time(time_step, limit_side);
    }
    
    if(set_time_result_ == TimeStatus::changed_forced)
        set_time_result_ = TimeStatus::changed;
    else
        set_time_result_ = TimeStatus::constant;
    
    // read all descriptors satisfying time.ge(input_time)
    update_history(time_step);
    check_initialized_region_fields_();

    //
    is_jump_time_=false;
    // set time_step on all regions
    // for regions that match type of the field domain
    for(const Region &reg: mesh()->region_db().get_region_set("ALL") ) {
    	auto rh = data_->region_history_[reg.idx()];

    	// Check regions with empty history, possibly set default.
    	if ( rh.empty()) continue;

        double last_time_in_history = rh.front().first;
        unsigned int history_size=rh.size();
        unsigned int i_history;
        ASSERT( time_step.ge(last_time_in_history) ).error("Setting field time back in history not fully supported yet!");

        // set history index
        if ( time_step.gt(last_time_in_history) ) {
            // in smooth time_step
            i_history=0;
        } else {
            // time_step .eq. input_time; i.e. jump time
            is_jump_time_=true;
            if (limit_side == LimitSide::right) {
                i_history=0;
            } else {
                i_history=1;
            }
        }
        i_history=min(i_history, history_size - 1);
        // This is a VERY old menaingless assert. Let see if we can live without it.
        // ASSERT(i_history >= 0).error("Empty field history.");

        // possibly update field pointer

        auto new_ptr = rh.at(i_history).second;
        if (new_ptr != region_fields_[reg.idx()]) {
            region_fields_[reg.idx()]=new_ptr;
            set_time_result_ = TimeStatus::changed;
        }
        // let FieldBase implementation set the time
        if ( new_ptr->set_time(time_step) )  set_time_result_ = TimeStatus::changed;

    }

    return changed();
}


template<int spacedim, class Value>
void Field<spacedim, Value>::copy_from(const FieldCommon & other) {
	ASSERT( flags().match(FieldFlag::equation_input))(other.name())(this->name())
	        .error("Can not copy to the non-input field.");

	// do not use copy if the field have its own input
	if ( flags().match(FieldFlag::declare_input)
	     && this->shared_->input_list_.size() != 0 ) return;

	if (typeid(other) == typeid(*this)) {
		auto  const &other_field = dynamic_cast<  Field<spacedim, Value> const &>(other);
		this->operator=(other_field);
	}
}



template<int spacedim, class Value>
void Field<spacedim, Value>::field_output(std::shared_ptr<OutputTime> stream, OutputTime::DiscreteSpace type)
{
	// currently we cannot output boundary fields
    ASSERT_LT( type, OutputTime::N_DISCRETE_SPACES ).error();
    this->compute_field_data( type, stream);
}



template<int spacedim, class Value>
FieldResult Field<spacedim,Value>::field_result( RegionSet region_set) const {

    FieldResult result_all = result_none;
    for(Region &reg : region_set) {
        auto f = region_fields_[reg.idx()];
        if (f) {
            FieldResult fr = f->field_result();
            if (result_all == result_none) // first region
                result_all = fr;
            else if (fr != result_all)
                result_all = result_other; // if results from individual regions are different
        } else return result_none; // if field is undefined on any region of the region set
    }

    if (result_all == result_constant && region_set.size() > 1)
        return result_other; // constant result for individual regions could be non-constant on the whole region set

    return result_all;

}


template<int spacedim, class Value>
std::string Field<spacedim,Value>::get_value_attribute() const
{
    int nrows = Value::NRows_;
    int ncols = Value::NCols_;
    string type = "Integer";
    if (std::is_floating_point<typename Value::element_type>::value)
        type = "Double";

    return fmt::format("{{ \"shape\": [ {}, {} ], \"type\": \"{}\", \"limit\": [ {}, {} ] }}",
    					nrows, ncols, type, this->limits().first, this->limits().second);
}


template<int spacedim, class Value>
void Field<spacedim,Value>::update_history(const TimeStep &time) {
	ASSERT_PTR( mesh() ).error("Null mesh pointer, set_mesh() has to be called before.\n");

    // read input up to given time
	double input_time;
    if (shared_->input_list_.size() != 0) {
        while( shared_->list_idx_ < shared_->input_list_.size()
        	   && time.ge( input_time = time.read_time( shared_->input_list_[shared_->list_idx_].find<Input::Tuple>("time") ) ) ) {

        	const Input::Record & actual_list_item = shared_->input_list_[shared_->list_idx_];
        	// get domain specification
        	RegionSet domain;
        	Input::Array domain_name_array;
        	unsigned int id;
			if (actual_list_item.opt_val("region", domain_name_array)) {
				std::vector<string> domain_names = mesh()->region_db().get_and_check_operands(domain_name_array);
				domain = mesh()->region_db().union_set(domain_names);

			} else if (actual_list_item.opt_val("rid", id)) {
				Region region;
				try {
					region = mesh()->region_db().find_id(id);
				} catch (RegionDB::ExcUniqueRegionId &e) {
					e << actual_list_item.ei_address();
					throw;
				}
				if (region.is_valid())
				    domain.push_back(region);
				else
				    THROW(RegionDB::ExcUnknownRegion() << RegionDB::EI_ID(id) );
			} else {
				THROW(ExcMissingDomain()
						<< actual_list_item.ei_address() );
			}

			// get field instance   
			for(auto rit = factories_.rbegin() ; rit != factories_.rend(); ++rit) {
				FieldBasePtr field_instance = (*rit)->create_field(actual_list_item, *this);
				if (field_instance)  // skip descriptors without related keys
				{
					// add to history
					ASSERT_EQ( field_instance->n_comp() , shared_->n_comp_);
					field_instance->set_mesh( mesh() );
					for(const Region &reg: domain) {
                        // if region history is empty, add new field
                        // or if region history is not empty and the input_time is higher, add new field
                        // otherwise (region history is not empty and the input_time is the same),
                        //      rewrite the region field
                        if( data_->region_history_[reg.idx()].size() == 0
                            || data_->region_history_[reg.idx()].back().first < input_time)
                        {
                            data_->region_history_[reg.idx()].push_front(
                                    HistoryPoint(input_time, field_instance));
                            //DebugOut() << "Update history" << print_var(this->name()) << print_var(reg.label()) << print_var(input_time);
                        }
                        else
                        {
                            data_->region_history_[reg.idx()].back() = 
                                    HistoryPoint(input_time, field_instance);
                        }
					}
					break;
				}
		    }

        	++shared_->list_idx_;
        }
    }
}

template<int spacedim, class Value>
void Field<spacedim,Value>::check_initialized_region_fields_() {
	ASSERT_PTR(mesh()).error("Null mesh pointer.");
    //if (shared_->is_fully_initialized_) return;

    // check there are no empty field pointers, collect regions to be initialized from default value
    RegionSet regions_to_init; // empty vector

    for(const Region &reg : mesh()->region_db().get_region_set("ALL") ) {
        RegionHistory &rh = data_->region_history_[reg.idx()];
        if ( rh.empty() ||	! rh[0].second)   // empty region history
        {
        		// test if check is turned on and control field is FieldConst
//                if (no_check_control_field_ && no_check_control_field_->is_constant(reg) ) {
//                	// get constant enum value
//                	auto elm = ElementAccessor<spacedim>(mesh(), reg);
//                	FieldEnum value = no_check_control_field_->value(elm.centre(),elm);
//                	// check that the value is in the disable list
//                    if ( std::find(shared_->no_check_values_.begin(), shared_->no_check_values_.end(), value)
//                             != shared_->no_check_values_.end() )
//                        continue;                  // the field is not needed on this region
//                }
            if (shared_->input_default_ != "") {    // try to use default
                regions_to_init.push_back( reg );
            } else {
                THROW( ExcMissingFieldValue() << EI_FieldInputName(input_name()) << EI_FieldName(name())
                        << EI_RegId(reg.id()) << EI_RegLabel(reg.label()) );
            }
        }
    }

    // possibly set from default value
    if ( regions_to_init.size() ) {
    	std::string region_list;
    	// has to deal with fact that reader can not deal with input consisting of simple values
    	string default_input=input_default();
    	auto input_type = get_input_type().make_instance().first;
        Input::ReaderToStorage reader( default_input, *input_type, Input::FileFormat::format_JSON );

        auto a_rec = reader.get_root_interface<Input::AbstractRecord>();
    	FieldAlgoBaseInitData init_data(input_name(), n_comp(), units(), limits(), flags());
        auto field_ptr = FieldBaseType::function_factory( a_rec , init_data );
        field_ptr->set_mesh( mesh() );
        for(const Region &reg: regions_to_init) {
    		data_->region_history_[reg.idx()]
    		                .push_front(HistoryPoint( 0.0, field_ptr) );
    		region_list+=" "+reg.label();
        }
        FieldCommon::messages_data_.push_back( MessageData(input_default(), name(), region_list) );

    }
    //shared_->is_fully_initialized_ = true;
}


template<int spacedim, class Value>
void Field<spacedim,Value>::add_factory(const std::shared_ptr<FactoryBase> factory) {
	factories_.push_back( factory );
}


template<int spacedim, class Value>
typename Field<spacedim,Value>::FieldBasePtr Field<spacedim,Value>::FactoryBase::create_field(Input::Record rec, const FieldCommon &field) {
	Input::AbstractRecord field_record;
	if (rec.opt_val(field.input_name(), field_record)) {
		FieldAlgoBaseInitData init_data(field.input_name(), field.n_comp(), field.units(), field.limits(), field.get_flags());
		return FieldBaseType::function_factory(field_record, init_data );
	}
	else
		return FieldBasePtr();
}


template<int spacedim, class Value>
bool Field<spacedim,Value>::FactoryBase::is_active_field_descriptor(const Input::Record &in_rec, const std::string &input_name) {
	return in_rec.find<Input::AbstractRecord>(input_name);
}




template<int spacedim, class Value>
void Field<spacedim,Value>::set_input_list(const Input::Array &list, const TimeGovernor &tg) {
    if (! flags().match(FieldFlag::declare_input)) return;

    // check that times forms ascending sequence
    double time,last_time=0.0;

    for (Input::Iterator<Input::Record> it = list.begin<Input::Record>();
					it != list.end();
					++it) {
    	for(auto rit = factories_.rbegin() ; rit != factories_.rend(); ++rit) {
			if ( (*rit)->is_active_field_descriptor( (*it), this->input_name() ) ) {
				shared_->input_list_.push_back( Input::Record( *it ) );
				time = tg.read_time( it->find<Input::Tuple>("time") );
				if (time < last_time) {
					THROW( ExcNonascendingTime()
							<< EI_Time(time)
							<< EI_Field(input_name())
							<< it->ei_address());
				}
				last_time = time;

				break;
			}
    	}
	}

}



template<int spacedim, class Value>
void Field<spacedim,Value>::set_output_data_cache(OutputTime::DiscreteSpace space_type, std::shared_ptr<OutputTime> stream) {
    typedef typename Value::element_type ElemType;

    auto output_cache_base = stream->prepare_compute_data<ElemType>(this->name(), space_type,
            (unsigned int)Value::NRows_, (unsigned int)Value::NCols_);
    output_data_cache_ = std::dynamic_pointer_cast<ElementDataCache<ElemType>>(output_cache_base);
}


template<int spacedim, class Value>
void Field<spacedim,Value>::compute_field_data(OutputTime::DiscreteSpace space_type, std::shared_ptr<OutputTime> stream) {
    std::shared_ptr<OutputMeshBase> output_mesh = stream->get_output_mesh_ptr();
    ASSERT(output_mesh);

    ASSERT_EQ(space_type, OutputTime::NATIVE_DATA);

    /* Copy data to array */
    std::shared_ptr< FieldFE<spacedim, Value> > field_fe_ptr = this->get_field_fe();

    if (field_fe_ptr) {
        auto native_output_data_base = stream->prepare_compute_data<double>(this->name(), space_type,
                (unsigned int)Value::NRows_, (unsigned int)Value::NCols_,
                typeid(field_fe_ptr->get_dofhandler()->ds()->fe()[0_d].get()).name(),  // should be used better solution of fe_type setting
                                                                                       // e.g. method 'name()' of FiniteElement and descendants
                field_fe_ptr->get_dofhandler()->max_elem_dofs());
        // try casting actual ElementDataCache
        auto native_output_data = std::dynamic_pointer_cast<ElementDataCache<double>>(native_output_data_base);
        field_fe_ptr->native_data_to_cache(*native_output_data);
    } else {
        WarningOut().fmt("Field '{}' of native data space type is not of type FieldFE. Output will be skipped.\n", this->name());
    }

    /* Set the last time */
    stream->update_time(this->time());

}


template<int spacedim, class Value>
void Field<spacedim,Value>::fill_data_value(const std::vector<int> &offsets)
{
    for (unsigned int i=0; i<offsets.size(); ++i) {
        if (offsets[i] == -1) continue; // skip empty value
        auto ret_value = Value::get_from_array(this->value_cache_, i);
        const Value &ele_value = Value( ret_value );
        output_data_cache_->store_value(offsets[i], ele_value.mem_ptr() );
    }
}


template<int spacedim, class Value>
void Field<spacedim,Value>::fill_observe_value(std::shared_ptr<ElementDataCacheBase> output_cache_base, const std::vector<int> &offsets)
{
    typedef typename Value::element_type ElemType;

    std::shared_ptr<ElementDataCache<ElemType>> observe_data_cache =
            std::dynamic_pointer_cast<ElementDataCache<ElemType>>(output_cache_base);

    for (unsigned int i=0; i<offsets.size(); ++i) {
        if (offsets[i] == -1) continue; // skip empty value
        auto ret_value = Value::get_from_array(this->value_cache_, i);
        const Value &ele_value = Value( ret_value );
        observe_data_cache->store_value(offsets[i], ele_value.mem_ptr() );
    }
}


template<int spacedim, class Value>
std::shared_ptr< FieldFE<spacedim, Value> > Field<spacedim,Value>::get_field_fe() {
	ASSERT_EQ(this->mesh()->region_db().size(), region_fields_.size()).error();
	//ASSERT(!this->shared_->bc_).error("FieldFE output of native data is supported only for bulk fields!");

	std::shared_ptr< FieldFE<spacedim, Value> > field_fe_ptr;

	bool is_fe = (region_fields_.size()>0); // indicate if FieldFE is defined on all bulk regions
	is_fe = is_fe && region_fields_[1] && (typeid(*region_fields_[1]) == typeid(FieldFE<spacedim, Value>));
	for (unsigned int i=3; i<2*this->mesh()->region_db().bulk_size(); i+=2)
		if (!region_fields_[i] || (region_fields_[i] != region_fields_[1])) {
			is_fe = false;
			break;
		}
	if (is_fe) {
		field_fe_ptr = std::dynamic_pointer_cast<  FieldFE<spacedim, Value> >( region_fields_[1] );
	}

	return field_fe_ptr;
}


template<int spacedim, class Value>
void Field<spacedim, Value>::cache_reallocate(const ElementCacheMap &cache_map, unsigned int region_idx) const {
    // Call cache_reinit of FieldAlgoBase descendant on appropriate region
	if (region_fields_[region_idx] != nullptr)
	    region_fields_[region_idx]->cache_reinit(cache_map);
}


template<int spacedim, class Value>
void Field<spacedim, Value>::cache_update(ElementCacheMap &cache_map, unsigned int region_patch_idx) const {
    unsigned int region_idx = cache_map.region_idx_from_chunk_position(region_patch_idx);
    if (region_fields_[region_idx] != nullptr) // skips bounadry regions for bulk fields and vice versa
        region_fields_[region_idx]->cache_update(value_cache_, cache_map, region_patch_idx);
}


template<int spacedim, class Value>
std::vector<const FieldCommon *> Field<spacedim, Value>::set_dependency(unsigned int i_reg) const {
   	if (region_fields_[i_reg] != nullptr) return region_fields_[i_reg]->set_dependency(*this->shared_->default_fieldset_);
   	else return std::vector<const FieldCommon *>();
}


template<int spacedim, class Value>
FieldValueCache<double> * Field<spacedim, Value>::value_cache() {
    return &value_cache_;
}


template<>
FieldValueCache<double> * Field<3, FieldValue<0>::Enum>::value_cache() {
    return nullptr;
}

template<>
FieldValueCache<double> * Field<3, FieldValue<0>::Integer>::value_cache() {
    return nullptr;
}


template<int spacedim, class Value>
const FieldValueCache<double> * Field<spacedim, Value>::value_cache() const {
    return &value_cache_;
}


template<>
const FieldValueCache<double> * Field<3, FieldValue<0>::Enum>::value_cache() const {
    return nullptr;
}

template<>
const FieldValueCache<double> * Field<3, FieldValue<0>::Integer>::value_cache() const {
    return nullptr;
}





#endif /* FIELD_IMPL_HH_ */
