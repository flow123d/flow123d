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

#include "fields/field_base.hh"
#include "io/output.h"

namespace IT=Input::Type;

/**
 *  Left and right time limit, used in the @p set_time() method.
 *  Assigned values allow to index an array.
 */
enum class LimitSide {
	left=0,
	right=1,
	unknown=2 	// undefined value
};


/**
 * @brief Common abstract parent of all Field<...> classes.
 *
 * We need common ancestor in order to keep a list of all fields in one EqData object and allow
 * collective operations like @p set_time or @p init_from_input.
 */
class FieldCommonBase {

public:
	TYPEDEF_ERR_INFO(EI_Time, double);
	TYPEDEF_ERR_INFO(EI_Field, std::string);
	DECLARE_INPUT_EXCEPTION(ExcNonascendingTime,
			<< "Non-ascending time: " << EI_Time::val << " for field " << EI_Field::qval << ".\n");
	DECLARE_INPUT_EXCEPTION(ExcMissingDomain,
			<< "Missing domain specification (region, r_id, or r_set) in fields descriptor:");
	DECLARE_EXCEPTION(ExcFieldMeshDifference,
			<< "Two copies of the field " << EI_Field::qval << "call set_mesh with different arguments.\n");



    /**
     *  Set name of the field, used for naming the field's key in EqData record.
     *  It can also be used to name a corresponding output data set, e.g. when outut the field into a VTK file.
     */
    FieldCommonBase &name(const string & name)
    { shared_->name_ = name; return *this;}
    /**
     * Set description of the field, used for description of corresponding key in documentation.
     */
    FieldCommonBase & desc(const string & desc)
    { shared_->desc_ = desc; return *this;}
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
    FieldCommonBase & input_default(const string &dflt)
    { shared_->default_ = dflt; return *this;}
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
    FieldCommonBase & units(const string & units)
    { shared_->units_ = units; return *this;}

    /**
     * For the fields returning "Enum", we have to pass the Input::Type::Selection object to
     * the field implementations.
     *
     * We must save raw pointer since selection may not be yet initialized (during static initialization phase).
     */
    FieldCommonBase & input_selection(const Input::Type::Selection *element_selection)
    {
      shared_->element_selection_=element_selection;
      return *this;
    }

    /**
     * Output data type used in the output() method. Can be different for different field copies.
     * one can choose between:
     * data constant on elements, linear data given in nodes, and discontinuous linear data.
     *
     * If not set explicitly by this method, the default value is Outputtime::ELEM_DATA
     */
    FieldCommonBase & output_type(OutputTime::RefType rt)
    { type_of_output_data_ = rt; return *this; }

    /**
     * Set number of components for run-time sized vectors. This is used latter when we construct
     * objects derived from FieldBase<...>.
     *
     * n_comp_ is constant zero for fixed values, this zero is set by Field<...> constructors
     */
    void n_comp( unsigned int n_comp)
    { shared_->n_comp_ = (shared_->n_comp_ ? n_comp : 0);}


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
     * with jump time. This method invalidate result of
     * @p changed() so it should be called just before @p set_time.
     * Can be different for different field copies.
     */
    void set_limit_side(LimitSide side)
    { limit_side_=side; }

    /**
     * Getters.
     */
    const std::string &name() const
    { return shared_->name_;}

    const std::string desc() const
    {return shared_->desc_;}

    const std::string &input_default() const
    { return shared_->default_;}

    const std::string &units() const
    { return shared_->units_;}

    OutputTime::RefType output_type() const
    { return type_of_output_data_; }

    bool is_bc() const
    { return shared_->bc_;}

    unsigned int n_comp() const
    { return shared_->n_comp_;}

    const Mesh * mesh() const
    { return shared_->mesh_;}

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
     * Returns input type for particular field instance, this is reference to a static member input_type of the corresponding @p FieldBase
     * class (i.e. with the same template parameters). This is used in FieldSet::make_field_descriptor_type.
     */
    virtual IT::AbstractRecord &get_input_type() =0;

    /**
     * Abstract method for initialization of the field on one region.
     */
    //virtual void set_from_input(const RegionSet &domain, const Input::AbstractRecord &rec) =0;

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
    virtual  bool set_time(const TimeGovernor &time=TimeGovernor()) =0;

    /**
     * Check that @p other is instance of the same Field<..> class and
     * perform assignment. Polymorphic copy.
     */
    virtual void copy_from(const FieldCommonBase & other) =0;

    /**
     * Output the field. Use type of output data given by @p type_of_output_data member.
     * The parameter @p output_rec is checked for key named by the field name. If the key exists its
     * string value is used to look for the OutputTime object of the same name, then the output of the field is performed.
     * If the key do not appear in the input, no output is done.
     */
    virtual void output(Input::Record output_rec) =0;


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
     * Virtual destructor.
     */
    virtual ~FieldCommonBase();


protected:
    /**
     * Private default constructor. Should be used only through
     * Field<...>
     */
    FieldCommonBase();

    /**
     * Private copy constructor. Should be used only through
     * Field<...>
     */
    FieldCommonBase(const FieldCommonBase & other);

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
	     * True for boundary fields.
	     */
	    bool bc_;
	    /**
	     * Number of components for fields that return variable size vectors. Zero in other cases.
	     */
	    unsigned int n_comp_;
	    /**
	     * Name of the particular field. Used to name the key in the Field list Record.
	     */
	    std::string name_;
	    /**
	     * Description of corresponding key in the Field list Record.
	     */
	    std::string desc_;
	    /**
	     * Units of the field values. Currently just a string description.
	     */
	    std::string units_;
	    /**
	     * For Enum valued fields this is the input type selection that should be used
	     * to read possible values of the field (e.g. for FieldConstant the key 'value' has this selection input type).
	     *
	     * Is empty selection for for non-enum values fields.
	     *
	     * In fact we must use raw pointer since selection may not be constructed yet (static variable).
	     */
	    const IT::Selection *element_selection_;
	    /**
	     * Possible default value of the field.
	     */
	    string default_;
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


	std::shared_ptr<SharedData> shared_;

    /**
     * Which value is returned for times where field is discontinuous.
     */
    LimitSide limit_side_;

    /**
     * Result of last set time method
     */
    enum class TimeStatus {
    	changed, 	//<  Field changed during last set time call.
    	constant,	//<  Field doesn't change.
    	unknown		//<  Before first call of set_time.
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
    OutputTime::RefType type_of_output_data_ = OutputTime::ELEM_DATA;

    /**
     * Maximum number of FieldBase objects we store per one region.
     */
    static const unsigned int history_length_limit_=3;


};







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
class Field : public FieldCommonBase {
public:

    typedef FieldBase<spacedim, Value> FieldBaseType;
    typedef std::shared_ptr< FieldBaseType > FieldBasePtr;
    typedef typename FieldBase<spacedim, Value>::Point Point;

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
    FieldBasePtr (*read_field_descriptor_hook)(Input::Record rec, const FieldCommonBase &field);

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
    static FieldBasePtr read_field_descriptor(Input::Record rec, const FieldCommonBase &field);

    /**
     * Check that whole field list is set, possibly use default values for unset regions
     * and call set_time for every field in the field list.
     *
     * Returns true if the field has been changed.
     */
    bool set_time(const TimeGovernor &time=TimeGovernor() ) override;

    /**
     * Check that other has same type and assign from it.
     */
    void copy_from(const FieldCommonBase & other) override;

    /**
     * Implementation of FieldCommonBase::output().
     */
    void output(Input::Record output_rec) override;


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
class MultiField : public FieldCommonBase {
public:
    //typedef FieldBase<spacedim, Value> SubFieldBaseType;
    typedef Field<spacedim, Value> SubFieldType;
    typedef Field<spacedim, typename FieldValue<spacedim>::Vector > TransposedField;

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
    void copy_from(const FieldCommonBase & other) override;

    /**
     * Implementation of @p FieldCommonBase::output().
     */
    void output(Input::Record output_rec) override;

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

    ASSERT(elm.region_idx().idx() < region_fields_.size(), "Region idx %u out of range %lu, field: %s\n",
           elm.region_idx().idx(), (unsigned long int) region_fields_.size(), name().c_str());
    ASSERT( region_fields_[elm.region_idx().idx()] ,
    		"Null field ptr on region id: %d, field: %s\n", elm.region().id(), name().c_str());

    region_fields_[elm.region_idx().idx()]->value_list(point_list,elm, value_list);
}







#endif /* FIELD_HH_ */
