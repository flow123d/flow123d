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
 * @file    field_common.cc
 * @brief   
 */

#include "fields/field_common.hh"

/****************************************************************************
 *  Implementation of FieldCommon
 */


FieldCommon::FieldCommon()
: shared_( std::make_shared<SharedData>() ),
  limit_side_(LimitSide::unknown),
  set_time_result_(TimeStatus::unknown),
  component_index_(std::numeric_limits<unsigned int>::max())
{
    shared_->bc_=false;
    shared_->input_default_="";
    shared_->n_comp_ = 0;
    shared_->mesh_ = nullptr;
    shared_->is_fully_initialized_=false;
}



FieldCommon::FieldCommon(const FieldCommon & other)
: shared_(other.shared_),
  limit_side_(LimitSide::unknown),
  set_time_result_(TimeStatus::unknown),
  component_index_(other.component_index_)
{
     flags_.add( FieldFlag::input_copy );
}



IT::Record FieldCommon::field_descriptor_record(const string& record_name) {
    return IT::Record(record_name, field_descriptor_record_decsription(record_name))
                     .declare_key("r_set", IT::String(), "Name of region set where to set fields.")
                     .declare_key("region", IT::String(), "Label of the region where to set fields. ")
                     .declare_key("rid", IT::Integer(0), "ID of the region where to set fields." )
                     .declare_key("time", IT::Double(0.0), IT::Default("0.0"),
                             "Apply field setting in this record after this time.\n"
                             "These times have to form an increasing sequence.")
					 .close();
}

const std::string FieldCommon::field_descriptor_record_decsription(const string& record_name) {
    return "Record to set fields of the equation.\n"
                "The fields are set only on the domain specified by one of the keys: 'region', 'rid', 'r_set'\n"
                "and after the time given by the key 'time'. The field setting can be overridden by\n"
                " any " + record_name + " record that comes later in the boundary data array.";
}



void FieldCommon::mark_input_times(TimeMark::Type mark_type) {
    if (! flags().match(FieldFlag::declare_input)) return;
    ASSERT_LESS( 0, shared_->input_list_.size());

    // pass through field descriptors containing key matching field name.
    double time;
    for( auto it = shared_->input_list_.begin();
         it != shared_->input_list_.end(); ++it) {
    	if (it->find<Input::Record>(input_name())) {
            time = it->val<double>("time"); // default time=0
            TimeGovernor::marks().add( TimeMark(time, mark_type | TimeGovernor::marks().type_input() ));
        }
    }
}



FieldCommon::~FieldCommon() {}
