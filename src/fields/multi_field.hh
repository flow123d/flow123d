/*
 * multi_field.hh
 *
 *  Created on: Feb 13, 2014
 *      Author: jb
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
 *  - implement set_complemented_vector_field
 *
 *  - problem with "input" methods, since Field works with AbstratRecord, the MultiField - However  - should use Array of AbstractRecords
 *    simplest solution - test that in EqDataBase and have more methods in FieldCommonBase, or somehow detach input handling from
 *    Fields
 *
 * Definition of MultiField must be in separate file.
 * In other case source file field.cc is too big and compiler can throw compile error.
 */
template<int spacedim, class Value>
class MultiField : public FieldCommon {
public:
    //typedef FieldBase<spacedim, Value> SubFieldBaseType;
    typedef Field<spacedim, Value> SubFieldType;
    typedef Field<spacedim, typename FieldValue<spacedim>::Vector > TransposedField;

    class MultiFieldFactory : public Field<spacedim, Value>::FactoryBase {
    public:
    	virtual typename Field<spacedim, Value>::FieldBasePtr create_field(Input::Record rec, const FieldCommon &field);
    };

    /**
     * Default constructor.
     */
    MultiField();

    /**
     * Returns input type of particular field instance, this is usually static member input_type of the corresponding FieldBase class (
     * with same template parameters), however, for fields returning "Enum" we have to create whole unique Input::Type hierarchy for
     * every instance since every such field use different Selection for initialization, even if all returns just unsigned int.
     */
    IT::AbstractRecord &get_input_type() override;

    IT::Record &get_multifield_input_type() override;

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
    bool set_time(const TimeGovernor &time) override;

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
    void output(OutputTime *stream) override;

    /**
     * Implementation of @p FieldCommonBase::is_constant().
     */
    bool is_constant(Region reg) override;

    /**
     * Virtual destructor.
     */
    inline virtual ~MultiField() {}

    /// Number of subfields that compose the multi-field.
    inline unsigned int size() const
    { return sub_fields_.size(); }

    /**
     * Allows set Field<dim, Vector> that can be used for alternative initialization in "transposed" form.
     */
    void set_complemented_vector_field( TransposedField &complemented);

    /**
     * Returns reference to the sub-field (component) of given index @p idx.
     */
    inline SubFieldType &operator[](unsigned int idx)
    {
    	ASSERT(idx < sub_fields_.size(), "Index of subfield is out of range.\n");
    	return sub_fields_[idx];
    }

private:
    std::vector< SubFieldType > sub_fields_;

    /// Helper class members, used only for input record
    SubFieldType sub_field_type_;
    TransposedField transposed_field_;
};


#endif /* MULTI_FIELD_HH_ */
