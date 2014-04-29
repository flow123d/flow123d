/*
 * field_set.hh
 *
 *  Created on: Mar 8, 2014
 *      Author: jb
 */

#ifndef FIELD_SET_HH_
#define FIELD_SET_HH_


#include <system/exceptions.hh>
#include <fields/field.hh>



/**
 * TODO: implementation robust against destroying fields before the FieldSet.
 */
class FieldSet {
public:
	DECLARE_EXCEPTION(ExcUnknownField, << "Field set has no field with name: " << FieldCommonBase::EI_Field::qval);

	/**
	 * Add an existing Field to the list. It stores just pointer to the field.
	 * Be careful to not destroy passed Field before the FieldSet.
	 *
	 * Using operator allows elegant setting and adding of a field to the field set:
	 * @code
	 * 		Field<...> init_quantity; // member of a FieldSet descendant
	 *
	 * 		field_set +=
	 * 			some_field
	 * 			.disable_where(type, {dirichlet, neumann}); // this must come first since it is not member of FieldCommonBase
	 * 			.name("init_temperature")
	 * 			.description("Initial temperature");
	 *
	 */
	FieldSet &operator +=(FieldCommonBase &add_field) {
		FieldCommonBase *found_field = field(add_field.name());
		if (found_field) {
			ASSERT(&add_field==found_field, "Another field of the same name exists when adding field: %s\n",
			        add_field.name().c_str());
		} else {
			field_list.push_back(&add_field);
			if (mesh_) add_field.set_mesh(*mesh_);
			if (!input_list_.is_empty()) add_field.set_input_list(input_list_);
			if (side_ != LimitSide::unknown) add_field.set_limit_side(side_);
		}
		return *this;
	}

	/**
	 * Add other FieldSet to current one.
	 */
	FieldSet &operator +=(const FieldSet &other) {
		for(auto field_ptr : other.field_list) this->operator +=(*field_ptr);
		return *this;
	}

	/**
	 * Make new FieldSet as a subset of *this. The new FieldSet contains fields with names given by the @p names parameter.
	 */
	FieldSet subset(std::vector<std::string> names) const {
		FieldSet set;
		for(auto name : names) set += (*this)[name];
		return set;
	}

	unsigned int size() const {
		return field_list.size();
	}

	/**
	 * Returns input type for a field descriptor, that can contain any of the fields in the set.
	 * Typical usage is from derived class, where we add fields in the constructor and make auxiliary temporary instance
	 * to get the record odf the field descriptor. Simplest example:
	 *
	 * @code
	 * class EqData : public FieldSet {
	 * public:
	 * 		// fields
	 * 		Field<..> field_a;
	 * 		Field<..> field_b
	 * 		EqData() {
	 * 			add(field_a);
	 * 			add(field_b);
	 * 		}
	 * }
	 *
	 * Input::Type::Record SomEquation::input_type=
	 * 		Record("SomeEquation","equation's description")
	 * 		.declare_key("data",Input::Type::Array(EqData().make_field_descriptor_type()),"List of field descriptors.");
	 * @endcode
	 */
    Input::Type::Record make_field_descriptor_type(const std::string &equation_name) const {
    	Input::Type::Record rec = FieldCommonBase::field_descriptor_record(equation_name + "_Data");
    	for(auto field : field_list) {
    		if (!field->is_just_copy()) {
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


    /**
     * Make Selection with strings for all field names in the FieldSet.
     */
    Input::Type::Selection make_output_field_selection(const string &name, const string &desc = "") {
    	namespace IT=Input::Type;
    	IT::Selection sel(name, desc);
    	int i=0;
    	// add value for each field excluding boundary fields
    	for( auto field : field_list)
    	{
    		if (!field->is_bc())
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
//    	sel.close();

    	return sel;
    }


    /**
     * Use @p FieldCommonBase::copy_from() to set field of the field set given by the first parameter @p dest_field_name.
     * The source field is given as the second parameter @p source. The field copies share same input descriptor list
     * and same instances of FieldBase classes but each copy can be set to different time and different limit side.
     */
    void set_field(const std::string &dest_field_name, FieldCommonBase &source) {
    	auto &field = (*this)[dest_field_name];
    	field.copy_from(source);
    	if (mesh_) ASSERT_EQUAL(mesh_, field.mesh() );
    	if (side_ != LimitSide::unknown) field.set_limit_side(side_);
    }



    /**
     * Return pointer to the field given by name @p field_name. Return nullptr if not found.
     */
    FieldCommonBase *field(const std::string &field_name) const {
        for(auto field : field_list)
            if (field->name() ==field_name) return field;
        return nullptr;
    }


    /**
     * Returns reference to the field given by @p field_name.
     * Throws if the field with given name is not found.
     */
    FieldCommonBase &operator[](const std::string &field_name) const {
        FieldCommonBase *found_field=field(field_name);
        if (found_field) return *found_field;

        THROW(ExcUnknownField() << FieldCommonBase::EI_Field(field_name));
        return *field_list[0]; // formal to prevent compiler warning
    }

    /**
     * Collective interface to @p FieldCommonBase::set_mesh().
     */
    void set_mesh(const Mesh &mesh) {
    	mesh_ = &mesh;
    	for(auto field : field_list) field->set_mesh(mesh);
    }

    /**
     * Collective interface to @p FieldCommonBase::set_mesh().
     */
    void set_input_list(Input::Array input_list) {
    	input_list_ = input_list;
    	for(auto field : field_list) field->set_input_list(input_list);
    }

    /**
     * Collective interface to @p FieldCommonBase::set_mesh().
     */
    void set_limit_side(LimitSide side) {
    	side_ = side;
    	for(auto field : field_list) field->set_limit_side(side);
    }
    /**
     * Collective interface to @p FieldCommonBase::mar_input_times().
     */
    void mark_input_times(TimeMark::Type mark_type) {
    	for(auto field : field_list) field->mark_input_times(mark_type);
    }
    /**
     * Collective interface to @p FieldCommonBase::set_mesh().
     */
    bool changed() const {
    	bool changed_all=false;
    	for(auto field : field_list) changed_all = changed_all || field->changed();
    	return changed_all;
    }

    /**
     * Collective interface to @p FieldCommonBase::set_mesh().
     */
    bool is_constant(Region reg) const {
    	bool const_all=true;
    	for(auto field : field_list) const_all = const_all && field->is_constant(reg);
    	return const_all;
    }

    /**
     * Collective interface to @p FieldCommonBase::set_mesh().
     */
    void set_time(const TimeGovernor &time) {
    	for(auto field : field_list) field->set_time(time);
    }

    /**
     * Collective interface to @p FieldCommonBase::output_type().
     * @param rt   Discrete function space (element, node or corner data).
     */
    void output_type(OutputTime::DiscreteSpace rt) {
    	for (auto field : field_list) field->output_type(rt);
	}

    /**
     * Collective interface to @p FieldCommonBase::output().
     */
    void output(OutputTime *stream) {
    	for(auto field : field_list) field->output(stream);
    }


    /**
     * Adds given field into list of fields for group operations on fields.
     * Parameters are: @p field pointer, @p name of the key in the input, @p desc - description of the key, and optional parameter
     * @p d_val with default value. This method is rather called through the macro ADD_FIELD
     */
    FieldCommonBase &add_field( FieldCommonBase *field, const string &name, const string &desc, const string & d_val="") {
    	*this += field->name(name).description(desc).input_default(d_val);
    	return *field;
    }

protected:


    /// List of all fields.
    std::vector<FieldCommonBase *> field_list;

    /// value set by last set_mesh(); set  the same to added fields
    const Mesh *mesh_ = nullptr;

    /// value set by last set_input_list(); set  the same to added fields
    Input::Array input_list_;

    /// value set by last set_time_limit(); set  the same to added fields
    LimitSide side_ = LimitSide::unknown;
 };


/**
 * (OBSOLETE)
 * Macro to simplify call of FieldSet::add_field method. Two forms are supported:
 *
 *
 *
 * ADD_FIELD(some_field, description);
 * ADD_FIELD(some_field, description, Default);
 *
 * The first form adds name "some_field" to the field member some_field, also adds description of the field. No default
 * value is specified, so the user must initialize the field on all regions (This is checked in the Field<..>::set_time method)
 *
 * The second form adds also default value to the field, that is Default(".."), or Default::read_time(), other default value specifications are
 * meaningless. The automatic conversion to FieldConst is used, e.g.  Default::("0.0") is automatically converted to
 * { TYPE="FieldConst", value=[ 0.0 ] } for a vector valued field, so you get zero vector on output on regions with default value.
 */

#define ADD_FIELD(name, ...)                   this->add_field(&name, string(#name), __VA_ARGS__)


#endif /* FIELD_SET_HH_ */
