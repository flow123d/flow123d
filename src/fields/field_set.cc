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
 * @file    field_set.cc
 * @brief   
 */

#include "fields/field_set.hh"
#include "system/sys_profiler.hh"



FieldSet &FieldSet::operator +=(FieldCommon &add_field) {
    FieldCommon *found_field = field(add_field.name());
    if (found_field) {
        ASSERT(&add_field==found_field, "Another field of the same name exists when adding field: %s\n",
                add_field.name().c_str());
    } else {
        field_list.push_back(&add_field);
    }
    return *this;
}



FieldSet &FieldSet::operator +=(const FieldSet &other) {
    for(auto field_ptr : other.field_list) this->operator +=(*field_ptr);
    return *this;
}



FieldSet FieldSet::subset(std::vector<std::string> names) const {
    FieldSet set;
    for(auto name : names) set += (*this)[name];
    return set;
}



FieldSet FieldSet::subset( FieldFlag::Flags::Mask mask) const {
    FieldSet set;
    for(auto field : field_list)
        if (field->flags().match(mask))  set += *field;
    return set;
}



Input::Type::Record FieldSet::make_field_descriptor_type(const std::string &equation_name) const {
    Input::Type::Record rec = Input::Type::Record(equation_name + "_Data",
    		FieldCommon::field_descriptor_record_description(equation_name + "_Data"))
    	.copy_keys(FieldCommon::field_descriptor_record(equation_name + "_Data_aux"));

    for(auto field : field_list) {
        if ( field->flags().match(FieldFlag::declare_input) ) {
            string description =  field->description() + " (($[" + field->units().format_latex() + "]$))";

            // Adding units is not so simple.
            // 1) It must be correct for Latex.
            // 2) It should be consistent with rest of documentation.
            // 3) Should be specified for all fields.
            //if (units != "") description+= " [" +field->units() + "]";

            // TODO: temporary solution, see FieldCommon::multifield_
            if (field->is_multifield()) {
            	rec.declare_key(field->input_name(), field->get_multifield_input_type(), description);
            } else {
            	rec.declare_key(field->input_name(), field->get_input_type(), description);
            }
        }

    }
    return rec.close();
}



Input::Type::Selection FieldSet::make_output_field_selection(const string &name, const string &desc)
{
    namespace IT=Input::Type;
    IT::Selection sel(name, desc);
    int i=0;
    // add value for each field excluding boundary fields
    for( auto field : field_list)
    {
        if ( !field->is_bc() && field->flags().match( FieldFlag::allow_output) )
        {
            string desc = "Output of the field " + field->name() + " (($[" + field->units().format_latex()+"]$))";
            if (field->description().length() > 0)
                desc += " (" + field->description() + ").";
            else
                desc += ".";
            sel.add_value(i, field->name(), desc);
            i++;
        }
    }

    return sel;
}



void FieldSet::set_field(const std::string &dest_field_name, FieldCommon &source)
{
    auto &field = (*this)[dest_field_name];
    field.copy_from(source);
}



FieldCommon *FieldSet::field(const std::string &field_name) const {
    for(auto field : field_list)
        if (field->name() ==field_name) return field;
    return nullptr;
}



FieldCommon &FieldSet::operator[](const std::string &field_name) const {
    FieldCommon *found_field=field(field_name);
    if (found_field) return *found_field;

    THROW(ExcUnknownField() << FieldCommon::EI_Field(field_name));
    return *field_list[0]; // formal to prevent compiler warning
}



bool FieldSet::changed() const {
    bool changed_all=false;
    for(auto field : field_list) changed_all = changed_all || field->changed();
    return changed_all;
}



bool FieldSet::is_constant(Region reg) const {
    bool const_all=true;
    for(auto field : field_list) const_all = const_all && field->is_constant(reg);
    return const_all;
}



void FieldSet::output(std::shared_ptr<OutputTime> stream) {
	START_TIMER("Fill OutputData");
    for(auto field : field_list)
        if ( !field->is_bc() && field->flags().match( FieldFlag::allow_output) )
            field->output(stream);
}



// OBSOLETE method
FieldCommon &FieldSet::add_field( FieldCommon *field, const string &name,
                                      const string &desc, const string & d_val) {
    *this += field->name(name).description(desc).input_default(d_val);
    return *field;
}



std::ostream &operator<<(std::ostream &stream, const FieldSet &set) {
    for(FieldCommon * field : set.field_list) {
        stream << *field
               << std::endl;
    }
    return stream;
}
