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

#include <boost/foreach.hpp>

#include "field.hh"
#include "field_algo_base.impl.hh"
#include "field_fe.hh"
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
: data_(std::make_shared<SharedData>())
{
	// n_comp is nonzero only for variable size vectors Vector, VectorEnum, ..
	// this invariant is kept also by n_comp setter
	shared_->n_comp_ = (Value::NRows_ ? 0 : 1);
	this->add_factory( std::make_shared<FactoryBase>() );

	this->multifield_ = false;
}


template<int spacedim, class Value>
Field<spacedim,Value>::Field(const string &name, bool bc)
: data_(std::make_shared<SharedData>())

{
		// n_comp is nonzero only for variable size vectors Vector, VectorEnum, ..
		// this invariant is kept also by n_comp setter
		shared_->n_comp_ = (Value::NRows_ ? 0 : 1);
		shared_->bc_=bc;
		this->name( name );
		this->add_factory( std::make_shared<FactoryBase>() );
		this->multifield_ = false;
}



template<int spacedim, class Value>
Field<spacedim,Value>::Field(unsigned int component_index, string input_name, string name, bool bc)
: data_(std::make_shared<SharedData>())
{
	// n_comp is nonzero only for variable size vectors Vector, VectorEnum, ..
	// this invariant is kept also by n_comp setter
	shared_->n_comp_ = (Value::NRows_ ? 0 : 1);
	this->set_component_index(component_index);
	this->name_ = (name=="") ? input_name : name;
	this->shared_->input_name_ = input_name;
    shared_->bc_ = bc;

	this->multifield_ = false;
}


template<int spacedim, class Value>
Field<spacedim,Value>::Field(const Field &other)
: FieldCommon(other),
  data_(other.data_),
  region_fields_(other.region_fields_),
  factories_(other.factories_)
{
	if (other.no_check_control_field_)
		no_check_control_field_ =  make_shared<ControlField>(*other.no_check_control_field_);

	this->multifield_ = false;
}


template<int spacedim, class Value>
Field<spacedim,Value> &Field<spacedim,Value>::operator=(const Field<spacedim,Value> &other)
{
	//OLD_ASSERT( flags().match( FieldFlag::input_copy )  , "Try to assign to non-copy field '%s' from the field '%s'.", this->name().c_str(), other.name().c_str());
	OLD_ASSERT(other.shared_->mesh_, "Must call set_mesh before assign to other field.\n");
	OLD_ASSERT( !shared_->mesh_ || (shared_->mesh_==other.shared_->mesh_),
	        "Assignment between fields with different meshes.\n");

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

	if (other.no_check_control_field_) {
		no_check_control_field_ =  make_shared<ControlField>(*other.no_check_control_field_);
	}

	return *this;
}



template<int spacedim, class Value>
it::Instance Field<spacedim,Value>::get_input_type() {
	return FieldBaseType::get_input_type_instance(shared_->input_element_selection_);
}



template<int spacedim, class Value>
it::Array Field<spacedim,Value>::get_multifield_input_type() {
	OLD_ASSERT(false, "This method can't be used for Field");

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
    OLD_ASSERT_LESS(reg.idx(), this->region_fields_.size());
    return this->region_fields_[reg.idx()];
}
*/


template <int spacedim, class Value>
bool Field<spacedim, Value>::is_constant(Region reg) {
	OLD_ASSERT(this->set_time_result_ != TimeStatus::unknown, "Unknown time status.\n");
	OLD_ASSERT_LESS(reg.idx(), this->region_fields_.size());
    FieldBasePtr region_field = this->region_fields_[reg.idx()];
    return (region_field && typeid(*region_field) == typeid(FieldConstant<spacedim, Value>));
}


template<int spacedim, class Value>
void Field<spacedim, Value>::set_field(
		const RegionSet &domain,
		FieldBasePtr field,
		double time)
{
	OLD_ASSERT(field, "Null field pointer.\n");

	OLD_ASSERT( mesh(), "Null mesh pointer, set_mesh() has to be called before set_field().\n");
    if (domain.size() == 0) return;

    OLD_ASSERT_EQUAL( field->n_comp() , n_comp());
    field->set_mesh( mesh() , is_bc() );

    HistoryPoint hp = HistoryPoint(time, field);
    for(const Region &reg: domain) {
    	RegionHistory &region_history = data_->region_history_[reg.idx()];
    	// insert hp into descending time sequence
    	OLD_ASSERT( region_history.size() == 0 || region_history[0].first < hp.first, "Can not insert smaller time %g then last time %g in field's history.\n",
    			hp.first, region_history[0].first );
    	region_history.push_front(hp);
    }
    set_history_changed();
}



template<int spacedim, class Value>
void Field<spacedim, Value>::set_field(
		const RegionSet &domain,
		const Input::AbstractRecord &a_rec,
		double time)
{
	FieldAlgoBaseInitData init_data(input_name(), n_comp(), units(), limits(), flags());
	set_field(domain, FieldBaseType::function_factory(a_rec, init_data), time);
}




template<int spacedim, class Value>
bool Field<spacedim, Value>::set_time(const TimeStep &time_step, LimitSide limit_side)
{
	OLD_ASSERT( mesh() , "NULL mesh pointer of field '%s'. set_mesh must be called before.\n",name().c_str());

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

    	// skip regions with no matching BC flag
    	if (reg.is_boundary() != is_bc()) continue;

    	// Check regions with empty history, possibly set default.
    	if ( rh.empty()) continue;

        double last_time_in_history = rh.front().first;
        unsigned int history_size=rh.size();
        unsigned int i_history;
        OLD_ASSERT(time_step.ge(last_time_in_history), "Setting field time back in history not fully supported yet!");

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
        OLD_ASSERT(i_history >= 0, "Empty field history.");
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
	ASSERT( flags().match(FieldFlag::equation_input))(other.name().c_str())(this->name().c_str())
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
void Field<spacedim, Value>::field_output(std::shared_ptr<OutputTime> stream)
{
	// currently we cannot output boundary fields
	if (!is_bc()) {
		unsigned int ids;
		const OutputTime::DiscreteSpace type = this->output_type();

		ASSERT_LT(type, OutputTime::N_DISCRETE_SPACES).error();

		OutputTime::DiscreteSpaceFlags flags = 1 << type;
	    for(ids=0; ids < OutputTime::N_DISCRETE_SPACES - 1; ids++)
	        if (flags & (1 << ids))
	        	this->compute_field_data( OutputTime::DiscreteSpace(ids), stream);
	    ids = OutputTime::NATIVE_DATA;
	    if (flags & (1 << ids))
	    	this->compute_native_data( OutputTime::DiscreteSpace(ids), stream);
	}
}


template<int spacedim, class Value>
void Field<spacedim, Value>::observe_output(std::shared_ptr<Observe> observe)
{
	typedef typename Value::element_type ElemType;

    if (observe->points().size() == 0) return;

    ElementDataCache<ElemType> &output_data = observe->prepare_compute_data<ElemType>(this->name(), this->time(),
    						(unsigned int)Value::NRows_, (unsigned int)Value::NCols_);

    unsigned int i_data=0;
    for(const ObservePoint &o_point : observe->points()) {
        unsigned int ele_index = o_point.element_idx();
        const Value &obs_value =
                        Value( const_cast<typename Value::return_type &>(
                                this->value(o_point.global_coords(),
                                        ElementAccessor<spacedim>(this->mesh(), ele_index, false)) ));
        ASSERT_EQ(output_data.n_elem(), obs_value.n_rows()*obs_value.n_cols()).error();
        output_data.store_value(i_data,  obs_value.mem_ptr());
        i_data++;
    }
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
	OLD_ASSERT( mesh(), "Null mesh pointer, set_mesh() has to be called before.\n");

    // read input up to given time
	double input_time;
    if (shared_->input_list_.size() != 0) {
        while( shared_->list_idx_ < shared_->input_list_.size()
        	   && time.ge( input_time = shared_->input_list_[shared_->list_idx_].val<double>("time") ) ) {

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
					OLD_ASSERT_EQUAL( field_instance->n_comp() , n_comp());
					field_instance->set_mesh( mesh() , is_bc() );
					for(const Region &reg: domain) {
						data_->region_history_[reg.idx()].push_front(
								HistoryPoint(input_time, field_instance)
						);
						//DebugOut() << "Update history" << print_var(this->name()) <<print_var(reg.label()) << print_var(input_time);
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
	OLD_ASSERT(mesh(), "Null mesh pointer.");
    //if (shared_->is_fully_initialized_) return;

    // check there are no empty field pointers, collect regions to be initialized from default value
    RegionSet regions_to_init; // empty vector

    for(const Region &reg : mesh()->region_db().get_region_set("ALL") )
        if (reg.is_boundary() == is_bc()) {      		// for regions that match type of the field domain
            RegionHistory &rh = data_->region_history_[reg.idx()];
        	if ( rh.empty() ||	! rh[0].second)   // empty region history
            {
        		// test if check is turned on and control field is FieldConst
                if (no_check_control_field_ && no_check_control_field_->is_constant(reg) ) {
                	// get constant enum value
                	auto elm = ElementAccessor<spacedim>(mesh(), reg);
                	FieldEnum value = no_check_control_field_->value(elm.centre(),elm);
                	// check that the value is in the disable list
                    if ( std::find(shared_->no_check_values_.begin(), shared_->no_check_values_.end(), value)
                             != shared_->no_check_values_.end() )
                        continue;                  // the field is not needed on this region
                }
                if (shared_->input_default_ != "") {    // try to use default
                    regions_to_init.push_back( reg );
                } else {
                	xprintf(UsrErr, "Missing value of the input field '%s' ('%s') on region ID: %d label: %s.\n",
                			input_name().c_str(), name().c_str(), reg.id(), reg.label().c_str() );
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
        field_ptr->set_mesh( mesh(), is_bc() );
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
void Field<spacedim,Value>::set_input_list(const Input::Array &list) {
    if (! flags().match(FieldFlag::declare_input)) return;

    // check that times forms ascending sequence
    double time,last_time=0.0;

    for (Input::Iterator<Input::Record> it = list.begin<Input::Record>();
					it != list.end();
					++it) {
    	for(auto rit = factories_.rbegin() ; rit != factories_.rend(); ++rit) {
			if ( (*rit)->is_active_field_descriptor( (*it), this->input_name() ) ) {
				shared_->input_list_.push_back( Input::Record( *it ) );
				time = it->val<double>("time");
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
void Field<spacedim,Value>::compute_field_data(OutputTime::DiscreteSpace space_type, std::shared_ptr<OutputTime> stream) {
	typedef typename Value::element_type ElemType;

	ASSERT(space_type != OutputTime::NATIVE_DATA).error();

    /* It's possible now to do output to the file only in the first process */
    if( stream->get_rank() != 0) {
        /* TODO: do something, when support for Parallel VTK is added */
        return;
    }

    ElementDataCache<ElemType> &output_data = stream->prepare_compute_data<ElemType>(this->name(), space_type,
    		(unsigned int)Value::NRows_, (unsigned int)Value::NCols_);


    /* Copy data to array */
    switch(space_type) {
    case OutputTime::NODE_DATA: {
        // set output data to zero
        vector<unsigned int> count(output_data.n_values(), 0);
        for(unsigned int idx=0; idx < output_data.n_values(); idx++)
            output_data.zero(idx);

        std::shared_ptr<OutputMesh> output_mesh = std::dynamic_pointer_cast<OutputMesh>( stream->get_output_mesh_ptr() );
        for(const auto & ele : *output_mesh )
        {
            std::vector<Space<3>::Point> vertices = ele.vertex_list();
            for(unsigned int i=0; i < ele.n_nodes(); i++)
            {
                unsigned int node_index = ele.node_index(i);
                const Value &node_value =
                        Value( const_cast<typename Value::return_type &>(
                                this->value(vertices[i],
                                            ElementAccessor<spacedim>(ele.orig_mesh(), ele.orig_element_idx(),false) ))
                             );
                output_data.add(node_index, node_value.mem_ptr() );
                count[node_index]++;
            }
        }

        // Compute mean values at nodes
        for(unsigned int idx=0; idx < output_data.n_values(); idx++)
            output_data.normalize(idx, count[idx]);
    }
    break;
    case OutputTime::CORNER_DATA: {
    	std::shared_ptr<OutputMeshDiscontinuous> output_mesh_disc
			= std::dynamic_pointer_cast<OutputMeshDiscontinuous>( stream->get_output_mesh_ptr(true) );
        for(const auto & ele : *output_mesh_disc )
        {
            std::vector<Space<3>::Point> vertices = ele.vertex_list();
            for(unsigned int i=0; i < ele.n_nodes(); i++)
            {
                unsigned int node_index = ele.node_index(i);
                const Value &node_value =
                        Value( const_cast<typename Value::return_type &>(
                        		this->value(vertices[i],
                                            ElementAccessor<spacedim>(ele.orig_mesh(), ele.orig_element_idx(),false) ))
                             );
                ASSERT_EQ(output_data.n_elem(), node_value.n_rows()*node_value.n_cols()).error();
                output_data.store_value(node_index, node_value.mem_ptr() );
            }
        }
    }
    break;
    case OutputTime::ELEM_DATA: {
    	std::shared_ptr<OutputMesh> output_mesh = std::dynamic_pointer_cast<OutputMesh>( stream->get_output_mesh_ptr() );
        for(const auto & ele : *output_mesh )
        {
            unsigned int ele_index = ele.idx();
            const Value &ele_value =
                        Value( const_cast<typename Value::return_type &>(
                        		this->value(ele.centre(),
                                            ElementAccessor<spacedim>(ele.orig_mesh(), ele.orig_element_idx(),false))
                                                                        )
                             );
            ASSERT_EQ(output_data.n_elem(), ele_value.n_rows()*ele_value.n_cols()).error();
            output_data.store_value(ele_index, ele_value.mem_ptr() );
        }
    }
    break;
    case OutputTime::NATIVE_DATA:
        //should not happen
    break;
    }

    /* Set the last time */
    stream->update_time(this->time());

}


template<int spacedim, class Value>
void Field<spacedim,Value>::compute_native_data(OutputTime::DiscreteSpace space_type, std::shared_ptr<OutputTime> stream) {
	ASSERT_EQ(space_type, OutputTime::NATIVE_DATA).error();

    /* It's possible now to do output to the file only in the first process */
    if( stream->get_rank() != 0) {
        /* TODO: do something, when support for Parallel VTK is added */
        return;
    }

    std::shared_ptr< FieldFE<spacedim, Value> > field_fe_ptr = this->get_field_fe();

    if (field_fe_ptr) {
        ElementDataCache<double> &output_data = stream->prepare_compute_data<double>(this->name(), space_type,
                (unsigned int)Value::NRows_, (unsigned int)Value::NCols_);
        field_fe_ptr->fill_data_to_cache(output_data);
    } else {
        WarningOut().fmt("Field '{}' of native data space type is not of type FieldFE. Output will be skipped.\n", this->name());
    }

    /* Set the last time */
    stream->update_time(this->time());

}


template<int spacedim, class Value>
std::shared_ptr< FieldFE<spacedim, Value> > Field<spacedim,Value>::get_field_fe() {
	ASSERT_EQ_DBG(this->mesh()->region_db().size(), region_fields_.size()).error();
	ASSERT(!this->shared_->bc_).error("FieldFE output of native data is supported only for bulk fields!");

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




#endif /* FIELD_IMPL_HH_ */
