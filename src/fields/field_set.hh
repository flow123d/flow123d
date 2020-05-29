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
 * @file    field_set.hh
 * @brief   
 */

#ifndef FIELD_SET_HH_
#define FIELD_SET_HH_


#include <iosfwd>                  // for ostream
#include <string>                  // for string
#include <vector>                  // for vector
#include "fields/field_common.hh"  // for FieldCommon, FieldCommon::EI_Field
#include "fields/field_flag.hh"    // for FieldFlag, FieldFlag::Flags
#include "fields/eval_subset.hh"   // for EvalSubset
#include "fields/eval_points.hh"   // for EvalPoints
#include "fields/field_value_cache.hh"
#include "input/accessors.hh"      // for Array
#include "input/type_record.hh"    // for Record
#include "io/output_time.hh"       // for OutputTime, OutputTime::DiscreteSpace
#include "system/exceptions.hh"    // for ExcStream, operator<<, DECLARE_EXC...
#include "system/flag_array.hh"    // for FlagArray<>::Mask
#include "tools/time_governor.hh"  // for TimeGovernor (ptr only), TimeStep
class Mesh;
class Region;



/**
 * @brief Container for various descendants of FieldCommonBase.
 *
 * Provides various collective operations.
 * Typical usage:
 *
 * class EqData : public FieldSet
 * {
 *      EqData() {
 *          *this += scalar_field
 *                  .name("scalar_field")
 *                  .description("Some description for input and output documentation.")
 *                  .input_default("{0.0}")
 *                  .units("m");
 *          *this += vector_field
 *                  .name("vector_field")
 *                  .description("Some description for input and output documentation.")
 *                  .units("m");
 *      }
 *
 *      Field<3, FieldValue<3>::Scalar> scalar_field;
 *      Field<3, FieldValue<3>::VectorFixed> vector_field;
 * };
 *
 * This way the fields are destructed just before their pointers stored in the FieldSet.
 *
 * TODO:
 * Some set_XY functions set also to the fields added to the FieldSet in future.
 * This behavior should be removed, since it is misleading in combination with mask subsets. If one set
 * something to mask subset, it does not influence fields added to the original field set even if
 * they match the mask of the subset.
 *
 */
class FieldSet : public FieldFlag {
public:
	DECLARE_EXCEPTION(ExcUnknownField, << "Field set has no field with name: " << FieldCommon::EI_Field::qval);

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
	 * 			.disable_where(type, {dirichlet, neumann}) // this must come first since it is not member of FieldCommonBase
	 * 			.name("init_temperature")
	 * 			.description("Initial temperature");
	 *
	 */
	FieldSet &operator +=(FieldCommon &add_field);

	/**
	 * Add other FieldSet to current one.
	 */
	FieldSet &operator +=(const FieldSet &other);

	/**
	 * Make new FieldSet as a subset of *this. The new FieldSet contains fields with names given by the @p names parameter.
	 */
	FieldSet subset(std::vector<std::string> names) const;

    /**
     * Make new FieldSet as a subset of *this.
     * The new FieldSet contains all fields that match given @p mask.
     */
    FieldSet subset( FieldFlag::Flags::Mask mask) const;

	/**
	 * Number of fields in the FieldSet.
	 */
	inline unsigned int size() const {
		return field_list.size();
	}

	/**
	 * Returns input type for a field descriptor, that can contain any of the fields in the set.
	 * Typical usage is from derived class, where we add fields in the constructor
	 * and make auxiliary temporary instance
	 * to get the record of the field descriptor.
	 * The returned Record has name equation_name + "_Data".
	 *
	 * Simplest example:
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
	 * 		.declare_key("data",Input::Type::Array(
	 * 		    EqData().make_field_descriptor_type("SomeEquation")),"List of field descriptors.");
	 * @endcode
	 *
	 */
    Input::Type::Record make_field_descriptor_type(const std::string &equation_name) const;

    /**
     * Make Selection with strings for all field names in the FieldSet.
     */
    //Input::Type::Selection make_output_field_selection(const string &name, const string &desc);

    /**
     * Use @p FieldCommonBase::copy_from() to set field of the field set given by the first parameter
     * @p dest_field_name. The source field is given as the second parameter @p source. The field
     * copies share the same input descriptor list and the same instances of FieldBase classes
     * but each copy can be set to different time and different limit side.
     *
     * See @p FieldCommonBase::copy_from documentation for details.
     */
    void set_field(const std::string &dest_field_name, FieldCommon &source);

    /**
     * Return pointer to the field given by name @p field_name. Return nullptr if not found.
     */
    FieldCommon *field(const std::string &field_name) const;

    /**
     * Returns reference to the field given by @p field_name.
     * Throws if the field with given name is not found.
     */
    FieldCommon &operator[](const std::string &field_name) const;

    /**
     * Collective interface to @p FieldCommonBase::set_components().
     * It is safe to call this for field sets containing also fields
     * with return value other then variable vector as long as all variable
     * vector fields should be set to the same number of components.
     */
    void set_components(const std::vector<string> &names) {
        for(FieldCommon *field : field_list) field->set_components(names);
    }
    /**
     * Collective interface to @p FieldCommonBase::set_mesh().
     */
    void set_mesh(const Mesh &mesh) {
        for(FieldCommon *field : field_list) field->set_mesh(mesh);
    }

    /**
     * Collective interface to @p FieldCommon::set_mesh().
     */
    void set_input_list(Input::Array input_list, const TimeGovernor &tg) {
        for(FieldCommon *field : field_list) field->set_input_list(input_list, tg);
    }

    /**
     * Collective interface to @p FieldCommonBase::flags_add().
     * @param mask   mask to set for all fields in the field set.
     */
    void flags_add( FieldFlag::Flags::Mask mask) {
        for(FieldCommon *field : field_list) field->flags_add(mask);
    }

    /**
     * Collective interface to @p FieldCommonBase::set_mesh().
     */
    bool set_time(const TimeStep &time, LimitSide limit_side);

    /**
     * Collective interface to @p FieldCommonBase::output_type().
     * @param rt   Discrete function space (element, node or corner data).
     */
    void output_type(OutputTime::DiscreteSpace rt) {
        for(FieldCommon *field : field_list) field->output_type(rt);
    }

    /**
     * Collective interface to @p FieldCommonBase::mark_input_times().
     */
    void mark_input_times(const TimeGovernor &tg) {
    	for(auto field : field_list) field->mark_input_times(tg);
    }

    /**
     * Collective interface to @p FieldCommonBase::set_mesh().
     */
    bool changed() const;

    /**
     * Collective interface to @p FieldCommonBase::set_mesh().
     */
    bool is_constant(Region reg) const;

    /**
     * Collective interface to @p FieldCommonBase::is_jump_time().
     */
    bool is_jump_time() const;

    /**
     * Collective interface to @p FieldCommon::cache_allocate().
     */
    void cache_allocate(std::shared_ptr<EvalPoints> eval_points) {
        for(auto field : field_list) field->cache_allocate(eval_points);
    }

    /**
     * Collective interface to @p FieldCommon::cache_update().
     */
    void cache_update(ElementCacheMap &cache_map) {
	    for(auto field : field_list) field->cache_update(cache_map);
    }

protected:

    /// List of all fields.
    std::vector<FieldCommon *> field_list;

    /**
     * Stream output operator
     */
    friend std::ostream &operator<<(std::ostream &stream, const FieldSet &set);
};



#endif /* FIELD_SET_HH_ */
