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
 * @file    field_common.hh
 * @brief   
 */

#ifndef FIELD_COMMON_HH_
#define FIELD_COMMON_HH_

#include <vector>
using namespace std;

#include "system/exceptions.hh"
#include "fields/field_values.hh"
#include "input/accessors.hh"
#include "input/type_generic.hh"
#include "tools/time_marks.hh"
#include "tools/time_governor.hh"

#include "fields/field_flag.hh"
#include "fields/unit_si.hh"
#include "io/output_time.hh"


class Region;

namespace IT=Input::Type;

/**
 *  Left and right time limit, used in the @p set_time() method.
 *  Assigned values allow to index an array.
 */
enum class LimitSide {
    left=0,
    right=1,
    unknown=2   // undefined value
};



/**
 * @brief Common abstract parent of all Field<...> classes.
 *
 * We need common ancestor in order to keep a list of all fields in one EqData object and allow
 * collective operations like @p set_time or @p init_from_input.
 */
class FieldCommon {

public:
    TYPEDEF_ERR_INFO(EI_Time, double);
    TYPEDEF_ERR_INFO(EI_Field, std::string);
    DECLARE_INPUT_EXCEPTION(ExcNonascendingTime,
            << "Non-ascending time: " << EI_Time::val << " for field " << EI_Field::qval << ".\n");
    DECLARE_INPUT_EXCEPTION(ExcMissingDomain,
            << "Missing domain specification (region, r_id, or r_set) in the field descriptor:");
    DECLARE_EXCEPTION(ExcFieldMeshDifference,
            << "Two copies of the field " << EI_Field::qval << "call set_mesh with different arguments.\n");



    /**
     *  Set name of the field. In fact there are two attributes set by this method.
     *
     *  The first is name used to identify the field as part of a FieldSet or MultiField objects.
     *  This name is permanent and can be set only by this method. Can be accessed by @p name() method.
     *  This name is also used at output.
     *
     *  The second is @p input_name_ that determines appropriate key name in the input field descriptor.
     *  This name is also set by this method, but is stored in the internal shared space which
     *  is overwritten during call of copy_from method or  assignment operator. Can be accessed by @p input_name() mathod.
     *
     */
    FieldCommon &name(const string & name)
    { name_=shared_->input_name_ = name;
      return *this;
    }
    /**
     * Set description of the field, used for description of corresponding key in documentation.
     */
    FieldCommon & description(const string & description)
    { shared_->input_description_ = description; return *this;}
    /**
     * Set default value for the field's key from which the default constant valued field will be constructed.
     *
     * During the first call of the @p set_time method, we check that the field is defined on all regions.
     * On regions where it is not set yet, we use given @p dflt string to get particular instance of
     * FieldBase<> (see @p check_initialized_region_fields_).
     * The default string is interpreted in the same way as if it appears in the input file
     * as the value of the field. In particular it can be whole record with @p TYPE of the field etc.
     * Most common choice is however mere constant.
     */
    FieldCommon & input_default(const string &input_default)
    { shared_->input_default_ = input_default; return *this;}
    /**
     * @brief Set basic units of the field.
     *
     * Currently, we use it only during output and we represents units just by a string.
     *
     * TODO:
     * Particular class for representing and conversion of various units would be more appropriate.
     * This can allow specification of the units on the inptu, automatic conversion and the same on the output.
     * Possibly this allow using Boost::Units library, however, it seems to introduce lot of boilerplate code.
     * But can increase correctness of the calculations.
     */
    FieldCommon & units(const UnitSI & units)
    { shared_->units_ = units; return *this;}

    /**
     * For the fields returning "Enum", we have to pass the Input::Type::Selection object to
     * the field implementations.
     *
     * We must save raw pointer since selection may not be yet initialized (during static initialization phase).
     */
    FieldCommon & input_selection(const Input::Type::Selection *element_selection)
    {
      shared_->input_element_selection_=element_selection;
      return *this;
    }

    /**
     * Output discrete space used in the output() method. Can be different for different field copies.
     * one can choose between:
     * data constant on elements, linear data given in nodes, and discontinuous linear data.
     *
     * If not set explicitly by this method, the default value is OutputTime::ELEM_DATA
     */
    FieldCommon & output_type(OutputTime::DiscreteSpace rt)
    { type_of_output_data_ = rt; return *this; }

    /**
     * Set given mask to the field flags, ignoring default setting.
     * Default setting is declare_input & equation_input & allow_output.
     */
    FieldCommon & flags(FieldFlag::Flags::Mask mask)
    { flags_ = FieldFlag::Flags(mask); return *this; }

    /**
     * Add given mask to the field flags.
     */
    FieldCommon & flags_add(FieldFlag::Flags::Mask mask)
    { flags().add(mask); return *this; }

    /**
     * Set vector of component names.
     * Set number of components for run-time sized vectors. This is used latter when we construct
     * objects derived from FieldBase<...>.
     *
     * n_comp_ is constant zero for fixed values, this zero is set by Field<...> constructors
     */
    void set_components(const std::vector<string> &names) {
        shared_->comp_names_ = names;
        shared_->n_comp_ = (shared_->n_comp_ ? names.size() : 0);
    }


    /**
     * Set internal mesh pointer.
     */
    virtual void set_mesh(const Mesh &mesh) {};
    /**
     * Set the data list from which field will read its input. It is list of "field descriptors".
     * When reading from the input list we consider only field descriptors containing key of
     * named by the field name. These field descriptors has to have times forming ascending sequence.
     *
     * The list is used by set_time method to set field on individual regions to actual FieldBase descendants.
     */
    void set_input_list(const Input::Array &list);

    /**
     * Set side of limit when calling @p set_time
     * with jump time, i.e. time where the field change implementation on some region.
     * Wee assume that implementations prescribe only smooth fields.
     * This method invalidate result of
     * @p changed() so it should be called just before @p set_time.
     * Can be different for different field copies.
     */
    void set_limit_side(LimitSide side)
    { this->limit_side_=side; }

    /**
     * Getters.
     */
    const std::string &input_name() const
    { return shared_->input_name_;}

    const std::string &name() const
    { return name_;}

    const std::string description() const
    {return shared_->input_description_;}

    const std::string &input_default() const
    { return shared_->input_default_;}

    const UnitSI &units() const
    { return shared_->units_;}

    OutputTime::DiscreteSpace output_type() const
    { return type_of_output_data_; }

    bool is_bc() const
    { return shared_->bc_;}

    unsigned int n_comp() const
    { return shared_->n_comp_;}

    const Mesh * mesh() const
    { return shared_->mesh_;}

    LimitSide limit_side() const
    { return limit_side_;}

    FieldFlag::Flags &flags()
    { return flags_; }

    /**
     * Returns time set by last call of set_time method.
     * Can be different for different field copies.
     */
    double time() const
    { return last_time_; }


    /**
     * Common part of the field descriptor. To get finished record
     * one has to add keys for individual fields. This is done automatically
     * using FieldSet::get_input_type().
     */
    static IT::Record field_descriptor_record(const string& record_name);

    /**
     * Create description of field descriptor record.
     */
    static const std::string field_descriptor_record_decsription(const string& record_name);

    /**
     * Returns input type for particular field instance, this is reference to a static member input_type of the corresponding @p FieldBase
     * class (i.e. with the same template parameters). This is used in FieldSet::make_field_descriptor_type.
     */
    virtual const IT::Instance &get_input_type() =0;

    /**
     * Returns input type for MultiField instance.
     * TODO: temporary solution, see @p multifield_
     */
    virtual IT::Record &get_multifield_input_type() =0;

    /**
     * Pass through the input array @p input_list_, collect all times where the field could change and
     * put appropriate time marks into global TimeMarks object.
     * Introduced time marks have both given @p mark_type and @p type_input() type.
     *
     * Further development:
     * - we have to distinguish "jump" times and "smooth" times
     */
    void mark_input_times(TimeMark::Type mark_type);

    /**
     * Abstract method to update field to the new time level.
     * Implemented by in class template Field<...>.
     *
     * Return true if the value of the field was changed on some region.
     * The returned value is also stored in @p changed_during_set_time data member.
     *
     * Default values helps when creating steady field. Note that default TimeGovernor constructor
     * set time to 0.0.
     *
     * Different field copies can be set to different times.
     */
    virtual  bool set_time(const TimeStep &time) =0;

    /**
     * Check that @p other is instance of the same Field<..> class and
     * perform assignment. Polymorphic copy.
     */
    virtual void copy_from(const FieldCommon & other) =0;

    /**
     * Output the field.
     * The parameter @p output_fields is checked for value named by the field name. If the key exists,
     * then the output of the field is performed. If the key do not appear in the input, no output is done.
     */
    virtual void output(std::shared_ptr<OutputTime> stream) =0;


    /**
     * If the field on given region @p reg exists and is of type FieldConstant<...> the method method returns true
     * otherwise it returns false.
     * Then one call ElementAccessor<spacedim>(mesh(), reg ) to construct an ElementAccessor @p elm
     * pointing to "virtual" element on which Field::value returns constant value.
     *
     * Current implementation use virtual functions and can be prohibitively slow if called for every element. If this
     * becomes necessary it is possible to incorporate such test into set_time method and in this method just return precomputed result.
     */
    virtual bool is_constant(Region reg) =0;

    /**
     * Returns true if set_time_result_ is not @p TimeStatus::constant.
     * Returns the same value as last set_time method.
     */
    bool changed() const
    {
        ASSERT( set_time_result_ != TimeStatus::unknown, "Invalid time status.\n");
        return ( (set_time_result_ == TimeStatus::changed) );
    }

    /**
     * Sets @p component_index_
     */
    void set_component_index(unsigned int idx)
    {
    	this->component_index_ = idx;
    }

    /**
     * Return @p multifield_ flag.
     * TODO: temporary solution
     */
    inline bool is_multifield() const
    {
    	return this->multifield_;
    }

    /**
     * Virtual destructor.
     */
    virtual ~FieldCommon();


protected:
    /**
     * Private default constructor. Should be used only through
     * Field<...>
     */
    FieldCommon();

    /**
     * Private copy constructor. Should be used only through
     * Field<...>
     */
    FieldCommon(const FieldCommon & other);

    /**
     * Invalidate last time in order to force set_time method
     * update region_fields_.
     */
    void set_history_changed()
    {
        last_time_ = -numeric_limits<double>::infinity();
    }

    /**
     * Setters for essential field properties.
     */
    /**
     *  Data shared among copies of the same field.
     *
     *  This allow field copies in different equations with different time setting, but
     *  sharing common input field descriptor array and common history.
     */
    struct SharedData {
    	/**
    	 * Empty constructor.
    	 */
    	SharedData() {};

        /**
         * True for boundary fields.
         */
        bool bc_;
        /**
         * Number of components for fields that return variable size vectors. Zero in other cases.
         */
        unsigned int n_comp_;
        /**
         * Names of field components.
         */
        std::vector< std::string > comp_names_;
        /**
         * Name of the particular field. Used to name the key in the Field list Record.
         */
        std::string input_name_;
        /**
         * Description of corresponding key in the Field list Record.
         */
        std::string input_description_;
        /**
         * Units of the field values. Currently just a string description.
         */
        UnitSI units_;
        /**
         * For Enum valued fields this is the input type selection that should be used
         * to read possible values of the field (e.g. for FieldConstant the key 'value' has this selection input type).
         *
         * Is empty selection for for non-enum values fields.
         *
         * In fact we must use raw pointer since selection may not be constructed yet (static variable).
         */
        const IT::Selection *input_element_selection_;
        /**
         * Possible default value of the field.
         */
        string input_default_;
        /**
         * Pointer to the mesh on which the field lives.
         */
        const Mesh *mesh_;

        /**
         * List of input field descriptors from which the field is set.
         */
        Input::Array input_list_;

        /**
         * Iterator to current input field descriptor.
         */
        Input::Iterator<Input::Record> list_it_;

        /**
         * True after check_initialized_region_fields_ is called. That happen at first call of the set_time method.
         */
        bool is_fully_initialized_;

        /**
         * For which values of an enum valued field we do not
         * check the field. User is responsible, that the value will not be called
         * on such regions.
         */
        std::vector<FieldEnum> no_check_values_;


    };

    /**
     * Name that identifies the field in the field_set. By default this is same as
     * shared_->input_name_.
     */
    std::string name_;

    /**
     * Data shared among copies of the same input field.
     */
    std::shared_ptr<SharedData> shared_;

    /**
     * Which value is returned for times where field is discontinuous.
     */
    LimitSide limit_side_;

    /**
     * Result of last set time method
     */
    enum class TimeStatus {
        changed,    //<  Field changed during last set time call.
        constant,   //<  Field doesn't change.
        unknown     //<  Before first call of set_time.
    };

    /// Status of @p history.
    TimeStatus set_time_result_;

    /**
     * Last set time. Can be different for different field copies.
     */
    double last_time_ = -numeric_limits<double>::infinity();

    /**
     * Output data type used in the output() method. Can be different for different field copies.
     */
    OutputTime::DiscreteSpace type_of_output_data_ = OutputTime::ELEM_DATA;

    /**
     * Specify if the field is part of a MultiField and which component it is
     */
    unsigned int component_index_;

    /**
     * Flag determining if object is Multifield or Field.
     * TODO: temporary solution, goal is to make these two classes to behave similarly
     */
    bool multifield_;

    /**
     * Maximum number of FieldBase objects we store per one region.
     */
    static const unsigned int history_length_limit_=3;

    /// Field flags. Default setting is "an equation input field, that can read from user input, and can be written to output"
    FieldFlag::Flags   flags_ = FieldFlag::declare_input & FieldFlag::equation_input & FieldFlag::allow_output;

    /**
     * Stream output operator
     */
    friend std::ostream &operator<<(std::ostream &stream, const FieldCommon &field) {

        string limit_side_str =
            (field.limit_side_ == LimitSide::left  ? "left" :
            (field.limit_side_ == LimitSide::right ? "right" :
              "unknown"));

        stream
        << "field name:" << field.name()
        << " limit side:" << limit_side_str
        << " n. comp.:" << field.n_comp()
        << " last time:" << field.last_time_;
        return stream;
    }
};







#endif /* FIELD_COMMON_HH_ */
