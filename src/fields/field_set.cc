/*
 * field_set.cc
 *
 *  Created on: Mar 8, 2014
 *      Author: jb
 */

#include "fields/field_set.hh"



FieldSet &FieldSet::operator +=(FieldCommonBase &add_field) {
    FieldCommonBase *found_field = field(add_field.name());
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
    Input::Type::Record rec = FieldCommonBase::field_descriptor_record(equation_name + "_Data");
    for(auto field : field_list) {
        if ( field->flags().match(FieldFlag::declare_input) ) {
            string units = field->units();
            string description =  field->description();

            // Adding units is not so simple.
            // 1) It must be correct for Latex.
            // 2) It should be consistent with rest of documentation.
            // 3) Should be specified for all fields.
            //if (units != "") description+= " [" +field->units() + "]";
            rec.declare_key(field->input_name(), field->get_input_type(), description);
        }

    }
    return rec;
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
            string desc = "Output of the field " + field->name(); //  + " [" + field->units() + "]";
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



void FieldSet::set_field(const std::string &dest_field_name, FieldCommonBase &source)
{
    auto &field = (*this)[dest_field_name];
    field.copy_from(source);
}



FieldCommonBase *FieldSet::field(const std::string &field_name) const {
    for(auto field : field_list)
        if (field->name() ==field_name) return field;
    return nullptr;
}



FieldCommonBase &FieldSet::operator[](const std::string &field_name) const {
    FieldCommonBase *found_field=field(field_name);
    if (found_field) return *found_field;

    THROW(ExcUnknownField() << FieldCommonBase::EI_Field(field_name));
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



void FieldSet::output(OutputTime *stream) {
    for(auto field : field_list)
        if ( !field->is_bc() && field->flags().match( FieldFlag::allow_output) )
            field->output(stream);
}



// OBSOLETE method
FieldCommonBase &FieldSet::add_field( FieldCommonBase *field, const string &name,
                                      const string &desc, const string & d_val) {
    *this += field->name(name).description(desc).input_default(d_val);
    return *field;
}



std::ostream &operator<<(std::ostream &stream, const FieldSet &set) {
    for(FieldCommonBase * field : set.field_list) {
        stream << *field
               << std::endl;
    }
    return stream;
}
