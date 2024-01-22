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

#include <iomanip>                  // for operator<<, setw, setfill, left
#include "fields/field_common.hh"
#include "input/accessors_impl.hh"  // for Record::val
#include "input/attribute_lib.hh"   // for Attribute
#include "input/storage.hh"         // for ExcStorageTypeMismatch
#include "tools/time_marks.hh"      // for TimeMark, TimeMark::Type, TimeMarks


/****************************************************************************
 *  Implementation of FieldCommon
 */


FieldCommon::FieldCommon()
: shared_( std::make_shared<SharedData>() ),
  set_time_result_(TimeStatus::unknown),
  is_jump_time_(true),
  component_index_(std::numeric_limits<unsigned int>::max())
{
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
                     .declare_key("time", TimeGovernor::get_input_time_type(), IT::Default("0.0"),
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
        time = tg.read_time( item.find<Input::Tuple>("time") ); // default time=0
        TimeGovernor::marks().add( TimeMark(time, mark_type | TimeGovernor::marks().type_input() ));
    }
}



bool FieldCommon::print_message_table(ostream& stream, std::string equation_name) {
	if (FieldCommon::messages_data_.size() == 0) return false;

	stream << endl << "Used default values of Fields for equation " << equation_name << ":" << endl;
	stream << std::setfill('-') << setw(100) << "" << endl;
	stream << std::setfill(' ') << " Field name" << setw(21) << "" << "Default value" << setw(7) << "" << "Apply on regions" << endl;
	for (std::vector<MessageData>::iterator it = FieldCommon::messages_data_.begin(); it < FieldCommon::messages_data_.end(); ++it) {
		stringstream ss;
		stream << " " << std::left << setw(30) << it->field_name_ << "" << " "
			   << setw(18) << it->default_value_ << "" << " " << it->region_list_ << endl;
	}
	stream << std::setfill('-') << setw(100) << "" << endl << endl;

	FieldCommon::messages_data_.clear();
	return true;
}



std::vector<FieldCommon::MessageData> FieldCommon::messages_data_ = std::vector<FieldCommon::MessageData>();



FieldCommon::~FieldCommon() {}
