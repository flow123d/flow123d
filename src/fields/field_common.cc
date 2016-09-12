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
  set_time_result_(TimeStatus::unknown),
  is_jump_time_(true),
  component_index_(std::numeric_limits<unsigned int>::max())
{
    shared_->bc_=false;
    shared_->input_default_="";
    shared_->n_comp_ = 0;
    shared_->mesh_ = nullptr;
    shared_->is_fully_initialized_=false;
}



FieldCommon::FieldCommon(const FieldCommon & other)
: name_(other.name_),
  shared_(other.shared_),
  set_time_result_(other.set_time_result_),
  last_time_(other.last_time_),
  last_limit_side_(other.last_limit_side_),
  is_jump_time_(other.is_jump_time_),
  component_index_(other.component_index_)
{
     flags_.add( FieldFlag::input_copy );
}



IT::Record FieldCommon::field_descriptor_record(const string& record_name) {
    return IT::Record(record_name, field_descriptor_record_description(record_name))
                     .declare_key("region", IT::Array( IT::String(), 1 ), "Labels of the regions where to set fields. ")
                     .declare_key("rid", IT::Integer(0), "ID of the region where to set fields.",
                             { {IT::Attribute::obsolete(),
                                     "\"Specification of the region by its ID is obsolete, will be removed in release 3.0.\\n"
                                     "Use region label declared in the Mesh record or default label 'region_<ID>'.\""} })
                     .declare_key("time", IT::Double(0.0), IT::Default("0.0"),
                             "Apply field setting in this record after this time.\n"
                             "These times have to form an increasing sequence.")
					 .close();
}

const std::string FieldCommon::field_descriptor_record_description(const string& record_name) {
    return "Record to set fields of the equation.\n"
                "The fields are set only on the domain specified by one of the keys: 'region', 'rid'\n"
                "and after the time given by the key 'time'. The field setting can be overridden by\n"
                " any " + record_name + " record that comes later in the boundary data array.";
}



void FieldCommon::mark_input_times(const TimeGovernor &tg) {
    if (! flags().match(FieldFlag::declare_input)) return;

    // pass through field descriptors containing key matching field name.
    TimeMark::Type mark_type = tg.equation_fixed_mark_type();
    double time;
    for( auto &item : shared_->input_list_) {
        time = item.val<double>("time"); // default time=0
        TimeGovernor::marks().add( TimeMark(time, mark_type | TimeGovernor::marks().type_input() ));
    }
}



FieldCommon::~FieldCommon() {}
