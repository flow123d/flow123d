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
#include "input/flow_attribute_lib.hh"
#include "fem/mapping_p1.hh"
#include "mesh/ref_element.hh"
#include "tools/bidirectional_map.hh"
#include "tools/unit_converter.hh"
#include <boost/algorithm/string/replace.hpp>
#include <queue>


FieldSet::FieldSet()
: mesh_(nullptr) {}


const Input::Type::Record & FieldSet::make_user_field_type(const std::string &equation_name) {
    static Field<3, FieldValue<3>::Scalar> scalar_field;
    static Field<3, FieldValue<3>::VectorFixed> vector_field;
    static Field<3, FieldValue<3>::TensorFixed> tensor_field;
    return Input::Type::Record( equation_name+":UserData", "Record to set fields of the equation: "+equation_name+".")
        .declare_key("name", Input::Type::String(), Input::Type::Default::obligatory(),
                     "Name of user defined field.")
        .declare_key("is_boundary", Input::Type::Bool(), Input::Type::Default("false"),
                     "Type of field: boundary or bulk.")
        .declare_key("scalar_field", scalar_field.get_input_type(),
                     "Instance of FieldAlgoBase ScalarField descendant.\n"
        		     "One of keys 'scalar_field', 'vector_field', 'tensor_field' must be set.\n"
        		     "If you set more than one of these keys, only first key is accepted.")
        .declare_key("vector_field", vector_field.get_input_type(),
                     "Instance of FieldAlgoBase VectorField descendant. See above for details.")
        .declare_key("tensor_field", tensor_field.get_input_type(),
                     "Instance of FieldAlgoBase TensorField descendant. See above for details.")
        .declare_key("unit", UnitConverter::get_input_type(), Input::Type::Default::optional(),
                     "Unit of the field values provided in the main input file, in the external file, or "
                     "by a function (FieldPython).")
		.close();
}

FieldSet &FieldSet::operator +=(FieldCommon &add_field) {
    FieldCommon *found_field = field(add_field.name());
    if (found_field) {
    	ASSERT_PERMANENT(&add_field==found_field)(add_field.name()).error("You cannot add field of the same name that exists in FieldSet!\n");
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
    set.set_mesh( *this->mesh_ );
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
    string rec_name = equation_name + ":Data";
    string desc = FieldCommon::field_descriptor_record_description(rec_name);
    Input::Type::Record rec = Input::Type::Record(rec_name, desc)
    	.copy_keys(FieldCommon::field_descriptor_record(rec_name));

    for(auto field : field_list) {
        if ( field->flags().match(FieldFlag::declare_input) ) {
            string description =  field->description() + " (($[" + field->units().format_latex() + "]$))";

            // Adding units is not so simple.
            // 1) It must be correct for Latex.
            // 2) It should be consistent with rest of documentation.
            // 3) Should be specified for all fields.
            //if (units != "") description+= " [" +field->units() + "]";

            // TODO: temporary solution, see FieldCommon::multifield_

            std::shared_ptr<Input::Type::TypeBase> field_type_ptr;
            if (field->is_multifield()) {
            	field_type_ptr = std::make_shared<Input::Type::Array>(field->get_multifield_input_type());
            } else {
                field_type_ptr = std::make_shared<Input::Type::Instance>(field->get_input_type());
            }
            ASSERT( field->units().is_def() )(field->input_name()).error("units not def.");
            Input::Type::TypeBase::attribute_map key_attributes = Input::Type::TypeBase::attribute_map(
                    { {FlowAttribute::field_unit(), field->units().json() },
                      {FlowAttribute::field_value_shape(), field->get_value_attribute()} }
           		);
            string default_val = field->input_default();
            if (default_val != "") {
            	boost::replace_all(default_val, "\"", "\\\"");
            	key_attributes[FlowAttribute::field_default_value()] = "\"" + default_val + "\"";
            }
            rec.declare_key(field->input_name(), field_type_ptr, Input::Type::Default::optional(), description, key_attributes);
        }

    }
    return rec.close();
}


/*
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
            DebugOut() << field->get_value_attribute();

            sel.add_value(i, field->name(), desc, { {FlowAttribute::field_value_shape(), field->get_value_attribute()} } );
            i++;
        }
    }

    return sel;
}
*/


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


bool FieldSet::set_time(const TimeStep &time, LimitSide limit_side) {
    bool changed_all=false;
    for(auto field : field_list) changed_all = field->set_time(time, limit_side) || changed_all;
    return changed_all;
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


bool FieldSet::is_jump_time() const {
    bool is_jump = false;
    for(auto field : field_list) is_jump = is_jump || field->is_jump_time();
    return is_jump;
}


void FieldSet::cache_update(ElementCacheMap &cache_map) {
    ASSERT_GT(region_field_update_order_.size(), 0).error("Variable 'region_dependency_list' is empty. Did you call 'set_dependency' method?\n");
    for (unsigned int i_reg_patch=0; i_reg_patch<cache_map.n_regions(); ++i_reg_patch) {
        for (const FieldCommon *field : region_field_update_order_[cache_map.region_idx_from_chunk_position(i_reg_patch)])
            field->cache_update(cache_map, i_reg_patch);
    }
}


void FieldSet::set_dependency(FieldSet &used_fieldset) {
    region_field_update_order_.clear();
    std::unordered_set<const FieldCommon *> used_fields;

    for (unsigned int i_reg=0; i_reg<mesh_->region_db().size(); ++i_reg) {
        for (FieldListAccessor f_acc : used_fieldset.fields_range()) {
            topological_sort( f_acc.field(), i_reg, used_fields );
        }
        used_fields.clear();
    }
}


void FieldSet::topological_sort(const FieldCommon *f, unsigned int i_reg, std::unordered_set<const FieldCommon *> &used_fields) {
    if (used_fields.find(f) != used_fields.end() ) return; // field processed
    used_fields.insert(f);
    auto dep_vec = f->set_dependency(i_reg); // vector of dependent fields
    for (auto f_dep : dep_vec) {
        topological_sort(f_dep, i_reg, used_fields);
    }
    region_field_update_order_[i_reg].push_back(f);
}


void FieldSet::add_coords_field() {
    *this += X_.name("X")
               .units(UnitSI().m())
               .input_default("0.0")
               .flags( FieldFlag::input_copy )
               .description("Coordinates field.");

    *this += depth_.name("d")
               .units(UnitSI().m())
               .input_default("0.0")
               .flags( FieldFlag::input_copy )
               .description("Depth field.");

    if (this->mesh_ != nullptr) {
        X_.set_mesh(*this->mesh_);
        depth_.set_mesh(*this->mesh_);
    }

    depth_.set_field_coords(&X_);
}


Range<FieldListAccessor> FieldSet::fields_range() const {
    auto bgn_it = make_iter<FieldListAccessor>( FieldListAccessor(field_list, 0) );
    auto end_it = make_iter<FieldListAccessor>( FieldListAccessor(field_list, field_list.size()) );
    return Range<FieldListAccessor>(bgn_it, end_it);
}


std::string FieldSet::print_dependency() const {
    ASSERT_GT(region_field_update_order_.size(), 0).error("Variable 'region_dependency_list' is empty. Did you call 'set_dependency' method?\n");
    std::stringstream s;
    for (auto reg_it : region_field_update_order_) {
        s << "\nregion_idx " << reg_it.first << ": ";
        for (auto f_it : reg_it.second) {
            s << f_it->name() << ", ";
        }
    }
    return s.str();
}


std::ostream &operator<<(std::ostream &stream, const FieldSet &set) {
    for(FieldCommon * field : set.field_list) {
        stream << *field
               << std::endl;
    }
    return stream;
}
