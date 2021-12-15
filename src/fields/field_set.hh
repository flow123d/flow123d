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
#include "fields/field.hh"
#include "fields/field_coords.hh"  // for FieldCoords
#include "fields/field_depth.hh"   // for FieldDepth
#include "fields/surface_depth.hh" // for SurfaceDepth
#include "mesh/range_wrapper.hh"
#include "tools/general_iterator.hh"
#include "input/accessors.hh"      // for Array
#include "input/type_record.hh"    // for Record
#include "io/output_time.hh"       // for OutputTime, OutputTime::DiscreteSpace
#include "system/exceptions.hh"    // for ExcStream, operator<<, DECLARE_EXC...
#include "system/flag_array.hh"    // for FlagArray<>::Mask
#include "tools/time_governor.hh"  // for TimeGovernor (ptr only), TimeStep
class Mesh;
class Region;
template <int spacedim, class Value> class FieldFormula;



/**
 * Accessor to vector of Fields holds in FieldSet.
 *
 * Class holds position to vector and allows iterate through all instances of Field class
 * and all components (SubFields) of MultiFields.
 *
 * Base methods:
 * - inc() - increment to next Field instance:
 *     Field - iterates to next item in field list
 *     MultiFields - iterates to next component of MultiField or if actual position is last component
 *                   jumps to next item in field list
 * - operator ->() - returns pointer to actual Field.
 */
class FieldListAccessor {
public:
    /// Default constructor
    FieldListAccessor()
    : field_idx_(0), field_component_idx_(0) {}

    /// Constructor
    FieldListAccessor(std::vector<FieldCommon *> field_list, unsigned int field_idx)
    : field_list_(field_list), field_idx_(field_idx), field_component_idx_(0) {}

    /// Iterates to next Field.
    inline void inc() {
        if (field_list_[field_idx_]->is_multifield()) {
            field_component_idx_++;
            if (field_component_idx_ == field_list_[field_idx_]->n_comp()) {
                field_idx_++;
                field_component_idx_ = 0;
            }
        } else {
        	field_idx_++;
        }
    }

    /// Getter for field_idx_
    inline unsigned int field_idx() const {
        return field_idx_;
    }

    /// Getter for field_component_idx_
    inline unsigned int field_component_idx() const {
        return field_component_idx_;
    }

	/// Returns pointer to actual field held by accessor
    FieldCommon * field() const {
        if (field_list_[field_idx_]->is_multifield())
            return field_list_[field_idx_]->get_component(field_component_idx_);
        else
            return field_list_[field_idx_];
    }

    /// Comparison of accessors.
	inline bool operator ==(const FieldListAccessor &other) {
		return this->field_idx_ == other.field_idx_ && field_component_idx_ == other.field_component_idx_;
	}

	inline bool operator !=(const FieldListAccessor &other) const {
		return this->field_idx_ != other.field_idx_ || field_component_idx_ != other.field_component_idx_;
	}

	/// Dereference operator simplify access to actual field held by accessor
    FieldCommon * operator ->() const {
        if (field_list_[field_idx_]->is_multifield())
            return field_list_[field_idx_]->get_component(field_component_idx_);
        else
            return field_list_[field_idx_];
    }

private:
    std::vector<FieldCommon *> field_list_;  ///< List of FieldCommon objects (combine Fields and MultiFields
    unsigned int field_idx_;                 ///< Index of actual Field in field_list
    unsigned int field_component_idx_;       ///< Index of subfield in MultiField (fo fields hold only value 0 that is not used)
};


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

	/// Default constructor.
	FieldSet();

	/**
    /**
     * @brief Declare input record type of field defined by user.
     */
    const Input::Type::Record & get_user_field(const std::string &equation_name);
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
    	this->mesh_ = &mesh;
        for(FieldCommon *field : field_list) field->set_mesh(mesh);
    }

    /**
     * Collective interface to @p FieldCommon::set_input_list().
     */
    void set_input_list(Input::Array input_list, const TimeGovernor &tg) {
        for(FieldCommon *field : field_list) field->set_input_list(input_list, tg);
    }

   /**
    * Fill data of user defined fields to user_fields_input_ map.
    */
   void set_user_fields_map(Input::Array input_list);

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
     * Collective interface to @p FieldCommon::recache_allocate().
     */
    void cache_reallocate(const ElementCacheMap &cache_map, FieldSet &used_fieldset) {
    	this->set_dependency(used_fieldset);
    	for (auto reg_it : region_field_update_order_) {
    	    unsigned int region_idx = reg_it.first;
    	    for (auto f_it : reg_it.second) {
    	        f_it->cache_reallocate(cache_map, region_idx);
    	    }
    	}
        //for(auto field : field_list) field->cache_reallocate(cache_map);
    }

    /**
     * Collective interface to @p FieldCommon::cache_update().
     */
    void cache_update(ElementCacheMap &cache_map);

    /**
     * Set reference of FieldSet to all instances of FieldFormula.
     */
    void set_dependency(FieldSet &used_fieldset);

    /**
     * Add coords field (X_) and depth field to field_list.
     *
     * We can't add this field automatically in constructor, because there is problem
     * in equation where we add one FieldSet to other.
     */
    void add_coords_field();

    /// Set surface depth object  to "d" field.
    inline void set_surface_depth(std::shared_ptr<SurfaceDepth> surface_depth) {
        depth_.set_surface_depth( surface_depth );
    }

    /// Returns range of Fields held in field_list
    Range<FieldListAccessor> fields_range() const;

    /// Returns pointer to mesh.
    inline const Mesh *mesh() const {
        return mesh_;
    }

    /// Return order of evaluated fields by dependency and region_idx.
    std::string print_dependency() const;

    /**
     * Return pointer to the field defined by user given by name @p field_name. Return nullptr if not found.
     */
    FieldCommon *user_field(const std::string &field_name, const TimeStep &time);


protected:

    /// Helper method sort used fields by dependency
    void topological_sort(const FieldCommon *f, unsigned int i_reg, std::unordered_set<const FieldCommon *> &used_fields);

    /// List of all fields.
    std::vector<FieldCommon *> field_list;

    /// List of fields defined by user.
    std::vector<FieldCommon *> user_field_list_;

    /// Pointer to the mesh.
    const Mesh *mesh_;

    /**
     * Holds vector of indices of fields in field_list sorted by dependency for every region.
     *
     * - first: index of region
     * - second: vector of indices of fields (corresponding to position in field_list vector)
     */
    std::map<unsigned int, std::vector<const FieldCommon *>> region_field_update_order_;

    // Default fields.
    // TODO derive from Field<>, make public, rename

    /// Field holds coordinates for computing of FieldFormulas
    FieldCoords X_;

    /// Field holds surface depth for computing of FieldFormulas
    FieldDepth depth_;

    /// Map assigns Input::Abstract to each field defined in optional Input::Array 'user_fields'
    std::unordered_map<std::string, Input::AbstractRecord> user_fields_input_;

    /**
     * Stream output operator
     */
    friend std::ostream &operator<<(std::ostream &stream, const FieldSet &set);

    template<int dim, class Val>
    friend class FieldFormula;
};



#endif /* FIELD_SET_HH_ */
