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
 * @file    multi_field.hh
 * @brief   
 */

#ifndef MULTI_FIELD_HH_
#define MULTI_FIELD_HH_

using namespace std;


#include "fields/field.hh"
#include "fields/field_common.hh"


namespace IT=Input::Type;

/**
 * @brief Class for representation of a vector of fields of the same physical quantity.
 *
 * When solving a system of same equations with the number of components given at runtime
 * (as in the case of transport equation for runtime given number of substances) we need means how to work with the whole
 * vector of fields at once. This is the aim of this class. It provides the interface given by the parent class @p FieldCommonBase,
 * but principally it is just a vector of Field<Value,dim> objects. The sub-fields or components of a @p MultiField are independent
 * objects, how ever the setters propagates the values from the MultiFields to the individual fields. The only exception is the
 * @p set_name method which in conjunction with @p MultiField::set_subfield_names can set unique name to each component.
 *
 * Template parameters are used for every subfield.
 *
 *  TODO:
 *  - general mechanism how to convert a Field< dim, Vector> to MultiField< dim, Value>
 *  - implement set_from_input
 *  - implement set_Time
 *
 *  - problem with "input" methods, since Field works with AbstratRecord, the MultiField - However  - should use Array of Abstracts
 *    simplest solution - test that in EqDataBase and have more methods in FieldCommonBase, or somehow detach input handling from
 *    Fields
 *
 * Definition of MultiField must be in separate file.
 * In other case source file field.cc is too big and compiler can throw compile error.
 */
template<int spacedim, class Value>
class MultiField : public FieldCommon {
public:
    typedef FieldAlgorithmBase<spacedim, Value> SubFieldBaseType;
    typedef Field<spacedim, Value> SubFieldType;
    typedef typename FieldAlgorithmBase<spacedim, Value>::Point Point;
    typedef FieldValue_<0,1,typename Value::element_type> MultiFieldValue;

    TYPEDEF_ERR_INFO( EI_MultiFieldName, const string );
    TYPEDEF_ERR_INFO( EI_Size, unsigned int );
    TYPEDEF_ERR_INFO( EI_ExpectedSize, unsigned int );
    DECLARE_INPUT_EXCEPTION( Exc_InvalidMultiFieldSize, << "Invalid size " << EI_Size::val
    		<< "of the MultiField " << EI_MultiFieldName::qval << ", expected size: " << EI_ExpectedSize::val );

    class MultiFieldFactory : public Field<spacedim, Value>::FactoryBase {
    public:
    	/// Constructor.
    	MultiFieldFactory(unsigned int index)
    	: index_(index) {}

    	virtual typename Field<spacedim, Value>::FieldBasePtr create_field(Input::Record rec, const FieldCommon &field);

    	bool is_active_field_descriptor(const Input::Record &in_rec, const std::string &input_name) override;

    	unsigned int index_;
    };

    /**
     * Default constructor.
     * bc indicates boundary multifield.
     */
    MultiField(bool bc = false);

    /**
     * Copy constructor.
     */
    MultiField(const MultiField &other);

    /**
     * Assignment operator. Same properties as copy constructor, but class member name_ is not copied.
     */
    MultiField &operator=(const MultiField &other);


    /**
     * Returns input type of particular field instance, this is usually static member input_type of the corresponding FieldBase class (
     * with same template parameters), however, for fields returning "Enum" we have to create whole unique Input::Type hierarchy for
     * every instance since every such field use different Selection for initialization, even if all returns just unsigned int.
     */
    IT::Instance get_input_type() override;

    IT::Array get_multifield_input_type() override;
    
    /**
     * By this method you can allow that the field need not to be set on regions (and times) where the given @p control_field is
     * FieldConstant and has value in given @p value_list. We check this in the set_time method. Through this mechanism we
     * can switch of e.g. boundary data fields according to the type of the boundary condition.
     */
    auto disable_where(
            const MultiField<spacedim, typename FieldValue<spacedim>::Enum > &control_field,
            const vector<FieldEnum> &value_list) -> MultiField &;

    /**
     * Abstract method to update field to the new time level.
     * Implemented by in class template Field<...>.
     *
     * Return true if the value of the field was changed on some region.
     * The returned value is also stored in @p changed_during_set_time data member.
     *
     * In first call initialize MultiField to the number of components given by the size of @p names
     * and use this vector  to name individual components. Should be called after the setters derived from
     * FieldCommonBase.
     */
    bool set_time(const TimeStep &time, LimitSide limit_side) override;

    /**
     * We have to override the @p set_mesh method in order to call set_mesh method for subfields.
     */
    void set_mesh(const Mesh &mesh) override;
    
    /**
     * Polymorphic copy. Check correct type, allows copy of MultiField or Field.
     */
    void copy_from(const FieldCommon & other) override;

    /**
     * Implementation of @p FieldCommonBase::output().
     */
    void output(std::shared_ptr<OutputTime> stream) override;

    /**
     * Implementation of FieldCommonBase::observe_output().
     */
    void observe_output(std::shared_ptr<Observe> observe) override;

    /**
     * Implementation of @p FieldCommonBase::is_constant().
     */
    bool is_constant(Region reg) override;

    /**
     * @brief Indicates special field states.
     *
     * Return possible values from the enum @p FieldResult, see description there.
     *
     * Only difference to Field<>::field_result is meaning of @p result_constant. Here we return this value if
     * all subfields are (possibly different) constants.
     *
     * Not used yet. Possibly tune behavior after usage.
     */
    FieldResult field_result( RegionSet region_set) const override;

    std::string get_value_attribute() const override;

    /**
     * Virtual destructor.
     */
    inline virtual ~MultiField() {}

    /// Number of subfields that compose the multi-field.
    inline unsigned int size() const
    { return sub_fields_.size(); }

    /**
     * Returns reference to the sub-field (component) of given index @p idx.
     */
    inline SubFieldType &operator[](unsigned int idx)
    {
    	OLD_ASSERT(idx < sub_fields_.size(), "Index of subfield in MultiField '%s' is out of range.\n", this->input_name().c_str());
    	return sub_fields_[idx];
    }
    
    /** 
     * Returns constant reference to the sub-field (component) of given index @p idx.
     */
    inline const SubFieldType &operator[](unsigned int idx) const
    {
    	OLD_ASSERT(idx < sub_fields_.size(), "Index of subfield in MultiField '%s' is out of range.\n", this->input_name().c_str());
        return sub_fields_[idx];
    }
    
    /**
     * Initialize components of MultiField.
     *
     * Must be call after setting components, mesh and limit side.
     */
    void setup_components();

    /**
     * Returns vector of value in one given point @p on an element given by ElementAccessor @p elm.
     * It returns reference to he actual value in order to avoid temporaries for vector and tensor values.
     */
    virtual typename MultiFieldValue::return_type value(const Point &p, const ElementAccessor<spacedim> &elm) const;

    /**
     * Returns std::vector of vector values in several points at once. The base class implements
     * trivial implementation using the @p value(,,) method. This is not optimal as it involves lot of virtual calls,
     * but this overhead can be negligible for more complex fields as Python of Formula.
     */
    virtual void value_list(const std::vector< Point >  &point_list, const  ElementAccessor<spacedim> &elm,
                             std::vector<typename MultiFieldValue::return_type>  &value_list) const;

    void set_input_list(const Input::Array &list) override;

private:
    /// Subfields (items) of MultiField
    std::vector< SubFieldType > sub_fields_;

    /// Full list of input field descriptors from which the subfields of MultiField are set.
    Input::Array full_input_list_;
    
    /**
     * If this pointer is set, turn off check of initialization in the
     * @p set_time method on the regions where the method @p get_constant_enum_value
     * of the control field returns value from @p no_check_values_. This
     * field is private copy, its set_time method is called from the
     * set_Time method of actual object.
     */
    const MultiField<spacedim, typename FieldValue<spacedim>::Enum > *no_check_control_field_;
};


#endif /* MULTI_FIELD_HH_ */
