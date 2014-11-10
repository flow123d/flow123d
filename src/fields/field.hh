/*
 * field.hh
 *
 *  Created on: Feb 13, 2014
 *      Author: jb
 */

#ifndef FIELD_HH_
#define FIELD_HH_

#include <memory>
using namespace std;

#include <boost/circular_buffer.hpp>

#include "system/exceptions.hh"
#include "input/accessors.hh"
#include "coupling/time_marks.hh"
#include "coupling/time_governor.hh"

#include "fields/field_common.hh"
#include "fields/field_algo_base.hh"
#include "fields/field_flag.hh"
#include "io/output_time.hh"


namespace IT=Input::Type;

/**
 * @brief Class template representing a field with values dependent on: point, element, and region.
 *
 * By "field" we mean a mapping of a a pair (Point, Time) to a @p Value, where
 * Point is from @p spacedim  dimensional ambient space, Time is real number (set by @p set_time method),
 * and @p Value type representing range of the field, which can be: real scalar, integer scalar (a discrete value),
 * real vector of fixed (compile time) size, real vector of runtime size, or a matrix of fixed dimensions.
 * Extensions to vectors or matrices of integers, or to variable tensors are possible. For vector and matrix values
 * we use classes provided by Armadillo library for linear algebra.
 *
 * This class assign particular fields (instances of descendants of FiledBase) to the regions. It keeps a table of pointers to fields for every possible bulk
 * region index (very same functionality, but for boundary regions is provided by @p BCField class). This class has interface very similar to  FiledBase, however
 * key methods @p value, and @p value_list are not virtual in this class by contrast these methods are inlined to minimize overhead for
 * simplest fields like FieldConstant.
 *
 * TODO: currently it works only for spacedim==3 since we have only mesh in 3D ambient space.
 *
 */
template<int spacedim, class Value>
class Field : public FieldCommon {
public:

    typedef FieldAlgorithmBase<spacedim, Value> FieldBaseType;
    typedef std::shared_ptr< FieldBaseType > FieldBasePtr;
    typedef typename FieldAlgorithmBase<spacedim, Value>::Point Point;

    static constexpr bool is_enum_valued = boost::is_same<typename Value::element_type, FieldEnum>::value;
    static const unsigned int space_dim = spacedim;


    /**
     * Pointer to function that creates an instance of FieldBase for
     * field with name @p field_name based on data in field descriptor @p rec.
     *
     * Default implementation in method @p read_field_descriptor  just reads key given by
     * @p field_name and creates instance using @p FieldBase<...>::function_factory.
     * Function should return empty SharedField (that is shared_ptr to FieldBase).
     *
     * Hooks are necessary to implement:
     * 1) backward compatibility with old BCD input files
     * 2) setting pressure values are piezometric head values
     */
    //FieldBasePtr (*read_field_descriptor_hook)(Input::Record rec, const FieldCommon &field);

    /**
     * Note for future:
     * We pass through parameter @p field information about field that holds the factory which are necessary
     * for interpreting user input and create particular field instance. It would be clearer to pass these information
     * when the factory is assigned to a field. Moreover some information may not be set to field at all but directly passed
     * to the factory.
     */
    class FactoryBase {
    public:
    	virtual FieldBasePtr create_field(Input::Record rec, const FieldCommon &field);
    };

    /**
     * Default constructor.
     *
     */
    Field();

    Field(const string &name, bool bc = false);

    /**
     * Copy constructor. Keeps shared history, declaration data, mesh.
     */
    Field(const Field &other);

    /**
     * Assignment operator. Same properties as copy constructor.
     *
     * Question: do we really need this, isn't copy constructor enough?
     */
    Field &operator=(const Field &other);


    /**
     * Returns reference to input type of particular field instance, this is static member @p input_type of the corresponding FieldBase class
     * (i.e. with same template parameters). However, for fields returning "Enum" we have to create whole unique Input::Type hierarchy using following method
     * @p meka_input_tree.
     * every instance since every such field use different Selection for initialization, even if all returns just unsigned int.
     */
    IT::AbstractRecord &get_input_type() override;


    /**
     * By this method you can allow that the field need not to be set on regions (and times) where the given @p control_field is
     * FieldConstant and has value in given @p value_list. We check this in the set_time method. Through this mechanism we
     * can switch of e.g. boundary data fields according to the type of the boundary condition.
     */
    auto disable_where(
    		const Field<spacedim, typename FieldValue<spacedim>::Enum > &control_field,
    		const vector<FieldEnum> &value_list) -> Field &;



    /**
     * Set mesh pointer and resize region arrays.
     *
     * Implements abstract method.
     */
    void set_mesh(const Mesh &mesh) override;


    /**
     * Direct read access to the table of Field pointers on regions.
     */
    //boost::shared_ptr< FieldBaseType > operator[] (Region reg);

    /**
     * Implementation of @p FieldCommonBase::is_constant().
     */
    bool is_constant(Region reg) override;


    /**
     * Assigns given @p field to all regions in given region set @p domain.
     * Field is added to the history with given time and possibly used in the next call of the set_time method.
     * Caller is responsible for correct construction of given field.
     *
     * Use this method only if necessary.
     *
     * Default time simplify setting steady fields.
     */
    void set_field(const RegionSet &domain, FieldBasePtr field, double time=0.0);

    /**
     * Same as before but the field is first created using FieldBase::function_factory(), from
     * given abstract record accessor @p a_rec.
     */
    void set_field(const RegionSet &domain, const Input::AbstractRecord &a_rec, double time=0.0);

    /**
     * Default implementation of @p read_field_descriptor_hook.
     *
     * Reads key given by @p field_name and creates the field instance using
     * @p FieldBase<...>::function_factory.
     */
    //static FieldBasePtr read_field_descriptor(Input::Record rec, const FieldCommon &field);

    void set_limit_side(LimitSide side) override
    { this->limit_side_=side; }

    /**
     * Check that whole field list is set, possibly use default values for unset regions
     * and call set_time for every field in the field list.
     *
     * Returns true if the field has been changed.
     */
    bool set_time(const TimeGovernor &time) override;

    /**
     * Check that other has same type and assign from it.
     */
    void copy_from(const FieldCommon & other) override;

    /**
     * Implementation of FieldCommonBase::output().
     */
    void output(OutputTime *stream) override;


    /**
     * Returns true, if field is currently set to a time in which it is discontinuous.
     */
    //bool is_jump_time();


    /**
     * Special field values spatially constant. Could allow optimization of tensor multiplication and
     * tensor or vector addition. field_result_ should be set in constructor and in set_time method of particular Field implementation.
     * We return value @p result_none, if the field is not initialized on the region of the given element accessor @p elm.
     */
    inline FieldResult field_result( ElementAccessor<spacedim> &elm) const;

    /**
     * Returns one value in one given point @p on an element given by ElementAccessor @p elm.
     * It returns reference to he actual value in order to avoid temporaries for vector and tensor values.
     */
    virtual typename Value::return_type const &value(const Point &p, const ElementAccessor<spacedim> &elm) const;

    /**
     * Returns std::vector of scalar values in several points at once. The base class implements
     * trivial implementation using the @p value(,,) method. This is not optimal as it involves lot of virtual calls,
     * but this overhead can be negligible for more complex fields as Python of Formula.
     */
    virtual void value_list(const std::vector< Point >  &point_list, const  ElementAccessor<spacedim> &elm,
                       std::vector<typename Value::return_type>  &value_list) const;

    /**
     * Set @p factory_base_ptr_
     */
    void set_factory_base_ptr(std::shared_ptr<FactoryBase> factory_base_ptr);

protected:
    /**
     * For fields returning "enum", i.e. with @p Value == FieldEnum, the input type (meaning whole input_Type tree of the field) depends on the
     * Input::Type::Selection object that represents particular C enum type. Therefore, we have to create whole tree for the selection
     * that was set through @p FieldBaseCommon::set_selection() method.
     */
    IT::AbstractRecord make_input_tree();




    /**
     * Read input into @p regions_history_ possibly pop some old values from the
     * history queue to keep its size less then @p history_length_limit_.
     */
    void update_history(const TimeGovernor &time);



    /**
     *  Check that whole field list (@p region_fields_) is set, possibly use default values for unset regions.
     */
    void check_initialized_region_fields_();

    /**************** Shared data **************/

    /// Pair: time, pointer to FieldBase instance
    typedef pair<double, FieldBasePtr> HistoryPoint;
    /// Nearest history of one region.
    typedef boost::circular_buffer<HistoryPoint> RegionHistory;

    struct SharedData {

        /**
         *  History for every region. Shared among copies.
         */
         std::vector< RegionHistory >  region_history_;
    };

    /**************** Data per copy **************/

    std::shared_ptr<SharedData> data_;

	/**
	 * If this pointer is set, turn off check of initialization in the
	 * @p set_time method on the regions where the method @p get_constant_enum_value
	 * of the control field returns value from @p no_check_values_. This
	 * field is private copy, its set_time method is called from the
	 * set_Time method of actual object.
	 */
    typedef Field<spacedim, typename FieldValue<spacedim>::Enum > ControlField;
	std::shared_ptr<ControlField>  no_check_control_field_;

    /**
     * Table with pointers to fields on individual regions.
     */
    std::vector< FieldBasePtr > region_fields_;

    std::shared_ptr<FactoryBase> factory_base_ptr_;



};





/**
 * Same as Field<...> but for boundary regions.
 */
template<int spacedim, class Value>
class BCField : public Field<spacedim, Value> {
public:
    BCField() : Field<spacedim,Value>("anonymous_bc", true) {}
};




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

    void set_limit_side(LimitSide side) override;

    /**
     * Abstract method to update field to the new time level.
     * Implemented by in class template Field<...>.
     *
     * Return true if the value of the field was changed on some region.
     * The returned value is also stored in @p changed_during_set_time data member.
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
     * Initialize MultiField to the number of components given by the size of @p names
     * and use this vector  to name individual components. Should be called after the setters derived from
     * FieldCommonBase.
     */
    void init( const vector<string> &names);

    /**
     * Allows set Field<dim, Vector> that can be used for alternative initialization in "transposed" form.
     */
    void set_complemented_vector_field( TransposedField &complemented);

    /**
     * Returns reference to the sub-field (component) of given index @p idx.
     */
    inline SubFieldType &operator[](unsigned int idx)
    { return sub_fields_[idx]; }

private:
    std::vector< SubFieldType > sub_fields_;
    std::vector< std::string > sub_names_;
};




/****************************************************************************************
 * Inlined methods of Field< ... >
 */

template<int spacedim, class Value>
inline typename Value::return_type const & Field<spacedim,Value>::value(const Point &p, const ElementAccessor<spacedim> &elm) const
{

    ASSERT(this->set_time_result_ != TimeStatus::unknown, "Unknown time status.\n");
    ASSERT(elm.region_idx().idx() < region_fields_.size(), "Region idx %u out of range %lu, field: %s\n",
           elm.region_idx().idx(), (unsigned long int) region_fields_.size(), name().c_str());
    ASSERT( region_fields_[elm.region_idx().idx()] ,
    		"Null field ptr on region id: %d, idx: %d, field: %s\n", elm.region().id(), elm.region_idx().idx(), name().c_str());
    return region_fields_[elm.region_idx().idx()]->value(p,elm);
}



template<int spacedim, class Value>
inline void Field<spacedim,Value>::value_list(const std::vector< Point >  &point_list, const ElementAccessor<spacedim> &elm,
                   std::vector<typename Value::return_type>  &value_list) const
{
    ASSERT(this->set_time_result_ != TimeStatus::unknown, "Unknown time status.\n");
    ASSERT(elm.region_idx().idx() < region_fields_.size(), "Region idx %u out of range %lu, field: %s\n",
           elm.region_idx().idx(), (unsigned long int) region_fields_.size(), name().c_str());
    ASSERT( region_fields_[elm.region_idx().idx()] ,
    		"Null field ptr on region id: %d, field: %s\n", elm.region().id(), name().c_str());

    region_fields_[elm.region_idx().idx()]->value_list(point_list,elm, value_list);
}







#endif /* FIELD_HH_ */
