/*
 * field_base.hh
 *
 *  Created on: Aug 31, 2012
 *      Author: jb
 */

/**
 * TODO:
 * - better tests:
 *   - common set of quantities with different kind of values (scalar, vector, tensor, discrete, ..),
 *     common points and elements for evaluation
 *   - for individual Field implementations have:
 *     - different input
 *     - possibly different EPETCT_EQ tests, but rather have majority common
 *
 *
 */


#ifndef FIELD_BASE_HH_
#define FIELD_BASE_HH_

#include <string>
#include <boost/type_traits.hpp>
#include <boost/shared_ptr.hpp>

#include "input/input_type.hh"
#include "input/accessors.hh"

#include "mesh/accessors.hh"
#include "mesh/point.hh"

#include "fields/field_values.hh"

template <int spacedim, class Value>
class FieldConstant;




/// Result type have sense only for larger Value types like vectors and tensors.
typedef enum  {
    result_none,    // field not set
    result_zero,    // zero scalar, vector, or tensor
    result_one,     // unit scalar (1.0), identity tensor
    result_other
} FieldResult;



/**
 * Base class for space-time function classes.
 */
template <int spacedim, class Value>
class FieldBase {
public:
       // expose template parameters
       typedef Value ValueType;
       static const unsigned int spacedim_=spacedim;


       /**
        * Kind of default constructor , with possible setting of the initial time.
        * Fields that returns variable size vectors accepts number of components @p n_comp.
        */
       FieldBase(unsigned int n_comp=0);

       /**
        * Returns template parameters as string in order to distinguish name of AbstractRecords
        * for initialization of different instances of the FieldBase template.
        */
       static std::string template_name();

       /**
        * Declaration of input type data member.
        */
       static Input::Type::AbstractRecord input_type;

       /**
        * Returns whole tree of input types for FieldBase with all descendants based on element input type (namely for FieldConstant)
        * given by element_input_type pointer. USE ONLY IF YOU CAN NOT USE
        * static member FieldBase<...>::input_type.
        */
       static Input::Type::AbstractRecord get_input_type(typename Value::ElementInputType *element_input_type=NULL);

       /**
        * This static method gets accessor to abstract record with function input,
        * dispatch to correct constructor and initialize appropriate function object from the input.
        * Returns shared pointer to  FunctionBase<>.
        */
       static boost::shared_ptr< FieldBase<spacedim, Value> >
           function_factory(const Input::AbstractRecord &rec, unsigned int n_comp=0);

       /**
        *  Function can provide way to initialize itself from the input data.
        */
       virtual void init_from_input(const Input::Record &rec);

       /**
        * Set new time value. Some Fields may and some may not implement time dependent values and
        * possibly various types of interpolation. There can not be unified approach to interpolation (at least not on this abstraction level)
        * since some fields (FieldFormula, FieldPython) provides naturally time dependent functions other fields like (FieldConstant, ...), however,
        * can be equipped by various time interpolation schemes. In future, we obviously need interpolation of variable order so that
        * we can use ODE integrators of higher order.
        *
        * The method returns true if the value of the field change in the new time step.
        */
       virtual bool set_time(double time);

       /**
        * Is used only by some Field implementations, but can be used to check validity of incoming ElementAccessor in value methods.
        *
        * Optional parameter @p boundary_domain can be used to specify, that the field will be evaluated only on the boundary part of the mesh.
        * TODO: make separate mesh for the boundary, then we can drop this parameter.
        */
       virtual void set_mesh(Mesh *mesh, bool boundary_domain);

       /**
        * Returns number of rows, i.e. number of components for variable size vectors. For values of fixed size returns zero.
        */
       unsigned int n_comp() const;

       /**
        * Special field values spatially constant. Could allow optimization of tensor multiplication and
        * tensor or vector addition. field_result_ should be set in constructor and in set_time method of particular Field implementation.
        */
        FieldResult field_result() const
        { return field_result_;}

       /**
        * Method for getting some information about next time where the function change its character.
        * Used to add appropriate TimeMarks.
        * TODO: think what kind of information we may need, is the next time value enough?
        */
       virtual double next_change_time()
       { ASSERT(0, "Not implemented yet."); return 0.0; }

       /**
        * Returns one value in one given point @p on an element given by ElementAccessor @p elm.
        * It returns reference to he actual value in order to avoid temporaries for vector and tensor values.
        *
        * This method just call the later one @p value(Point, ElementAccessor, Value) and drops the FieldResult.
        *
        * Usual implementation of this method fills @p member r_value_ through unified envelope @p value_ as general tensor
        * and then returns member @p r_value_. However, some particular Fields may have result stored elsewhere, in such a case
        * the reference to the result can be returned directly without using the member @p value_. Keeping such wide space for optimization
        * has drawback in slow generic implementation of the @p value_list method that fills whole vector of values for vector of points.
        * Its generic implementation has to copy all values instead of directly store them into the vector of result values.
        *
        * So the best practice when implementing @p value and @value_list methods in particular FieldBase descendant is
        * implement some thing like value(point, elm, Value::return_type &value) and using
        *  s having in part
        *
        */
       virtual typename Value::return_type const &value(const Point<spacedim> &p, const ElementAccessor<spacedim> &elm)=0;

       /**
        * Pure virtual method. At least this has to be implemented by descendants.
        * Returns one value in one given point. ResultType can be used to avoid some costly calculation if the result is trivial.
        */
       //virtual FieldResult value(const Point<spacedim> &p, ElementAccessor<spacedim> &elm, typename Value::return_type &value) =0;

       /**
        * Returns std::vector of scalar values in several points at once. The base class implements
        * trivial implementation using the @p value(,,) method. This is not optimal as it involves lot of virtual calls,
        * but this overhead can be negligible for more complex fields as Python of Formula.
        */
       virtual void value_list(const std::vector< Point<spacedim> >  &point_list, const ElementAccessor<spacedim> &elm,
                          std::vector<typename Value::return_type>  &value_list)=0;

       /**
        * Virtual destructor.
        */
       virtual ~FieldBase() {}


protected:
       /// Actual time level; initial value is -infinity.
       double time_;
       /// Last value, prevents passing large values (vectors) by value.
       Value value_;
       typename Value::return_type r_value_;
       /// Indicator of particular values (zero, one) constant over space.
       FieldResult field_result_;
};




/**
 * @brief Common abstract parent of all Field<...> classes.
 *
 * We need common ancestor in order to keep a list of all fields in one EqData object and allow
 * collective operations like @p set_time or @p init_from_input.
 */
class FieldCommonBase {
public:
    /**
     * Constructor, we denote if this is bulk or bc field.
     */
    FieldCommonBase(bool bc);
    /**
     *  Set name of the field, used for naming the field's key in EqData record.
     *  It can also be used to name a corresponding output data set, e.g. when outut the field into a VTK file.
     */
    inline void set_name(const string & name)   { name_ = name;}
    /**
     * Set description of the field, used for description of corresponding key in documentation.
     */
    inline void set_desc(const string & desc)   { desc_ = desc; }
    /**
     * Set default value for the field's key from which the default constant valued field will be constructed.
     * During the first call of the @p set_time method, the table that assign particular fields to the regions
     * is checked for NULL pointers that corresponds to the regions without explicit field specification from the input.
     * On these regions we use the default constant field.
     */
    inline void set_default(const IT::Default &dflt)    { default_ = dflt; }
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
    inline void set_units(const string & units)         { units_ = units; }
    /**
     * Set number of components for run-time sized vectors. This is used latter when we construct
     * objects derived from FieldBase<...>.
     */
    inline void set_n_comp( unsigned int n_comp)        { n_comp_ = n_comp; }
    /**
     * For the fields returning "Enum", we have to pass the Input::Type::Selection object to
     * the field implementations.
     */
    inline void set_selection( Input::Type::Selection *element_selection)   { element_selection_=element_selection;}
    /**
     * Set internal mesh pointer.
     */
    virtual inline void set_mesh(Mesh *mesh)                    { mesh_=mesh; }

    /**
     * Getters.
     */
    const std::string &name() const;
    const std::string desc() const;
    const IT::Default &get_default() const;
    const std::string &units() const;
    double time() const;
    bool is_bc() const;
    bool is_enum_valued() const;
    unsigned int n_comp() const;
    Mesh * mesh() const;
    bool changed() const;

    /**
     * Returns input type for particular field instance, this is reference to a static member input_type of the corresponding @p FieldBase
     * class (i.e. with the same template parameters). This is used by EqDataBase::generic_input_type to construct the Input::Type::Record for
     * the bc_data_list or bulk_data_list.
     */
    virtual IT::AbstractRecord &get_input_type() =0;
    /**
     * Returns whole input type tree for FieldBase class returning "Enum".  We have to create whole unique Input::Type hierarchy for
     * every instance since every such field use different Selection for initialization, even if all returns just unsigned int.
     * The Input::Type::Selection object has to be set by the  @p set_selection method since @p make_input_tree is called by
     * @p EqDataBase::generic_input_type, where the information about  the Selection is not available.
     */
    virtual IT::AbstractRecord make_input_tree() =0;

    /**
     * Abstract method for initialization of the field on one region.
     */
    virtual void set_from_input(const RegionSet &domain, const Input::AbstractRecord &rec) =0;

    /**
     * Abstract method to update field to the new time level.
     * Implemented by in class template Field<...>.
     *
     * Return true if the value of the field was changed on some region.
     * The returned value is also stored in @p changed_during_set_time data member.
     */
    virtual bool set_time(double time) =0;

    /**
     * Virtual destructor.
     */
    virtual ~FieldCommonBase();

    /**
     * Is true if the values of the field has changed during last set_time() call.
     */
    bool changed_during_set_time;

protected:
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
     * True for boundary fields.
     */
    bool bc_;
    /**
     * Number of components for fields that return variable size vectors. Zero in other cases.
     */
    unsigned int n_comp_;
    /**
     * For Enum valued fields this points to the input type selection that should be used
     * to read possible values of the field (e.g. for FieldConstant the key 'value' has this selection input type).
     */
    IT::Selection *element_selection_;
    /**
     * Possible default value of the field. This implies TYPE=FieldConstant.
     */
    IT::Default default_;
    /**
     * Is true if the value returned by the field is based on Enum
     *  (i.e. constant value is initialized by some Input::Type::Selection)
     */
    bool enum_valued_;
    /**
     * Pointer to the mesh on which the field lives.
     */
    Mesh *mesh_;

    /**
     * Set by other methods (namely set_field() and set_from_input()) that modify the field before the set_time is called.
     */
    bool changed_from_last_set_time_;

    /// Time of the last set time.
    double last_set_time_;

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
 *
 */
template<int spacedim, class Value>
class Field : public FieldCommonBase {
public:
    typedef FieldBase<spacedim, Value> FieldBaseType;

    /**
     * Default constructor.
     *
     */
    Field();

    /**
     * Returns reference to input type of particular field instance, this is static member @p input_type of the corresponding FieldBase class
     * (i.e. with same template parameters). However, for fields returning "Enum" we have to create whole unique Input::Type hierarchy using following method
     * @p meka_input_tree.
     * every instance since every such field use different Selection for initialization, even if all returns just unsigned int.
     */
    IT::AbstractRecord &get_input_type();

    /**
     * For fields returning "enum", i.e. with @p Value == FieldEnum, the input type (meaning whole input_Type tree of the field) depends on the
     * Input::Type::Selection object that represents particular C enum type. Therefore, we have to create whole tree for the selection
     * that was set through @p FieldBaseCommon::set_selection() method.
     */
    IT::AbstractRecord make_input_tree();

    /**
     * By this method you can allow that the field need not to be set on regions (and times) where the given @p control_field is
     * FieldConstant and has value in given @p value_list. We check this in the set_time method. Through this mechanism we
     * can switch of e.g. boundary data fields according to the type of the boundary condition.
     */
    void disable_where(const Field<spacedim, typename FieldValue<spacedim>::Enum > *control_field, const vector<FieldEnum> &value_list);

    /**
     * Direct read access to the table of Field pointers on regions.
     */
    boost::shared_ptr< FieldBaseType > operator[] (Region reg);

    /**
     * If the field on given region @p reg exists and is of type FieldConstant<...> the method method returns true and sets
     * given ElementAccessor @p elm to "virtual" element on which Field::value returns constant value.
     * Otherwise it returns false and invalidate @p elm (ElementAccessor::is_valid() returns false).
     */
    bool get_const_accessor(Region reg, ElementAccessor<spacedim> &elm);

    /**
     * Initialize field of region @p reg from input accessor @p rec. At first usage it allocates
     * table of fields according to the @p bulk_size of the RegionDB. RegionDB is automatically closed.
     */
    void set_from_input(const RegionSet &domain, const Input::AbstractRecord &rec);

    /**
     * Assigns @p field to the given region @p reg. Caller is responsible for correct construction of given field.
     * Use this method only if necessary.
     */
    void set_field(const RegionSet &domain, boost::shared_ptr< FieldBaseType > field);

    /**
     * Check that whole field list is set, possibly use default values for unset regions
     *  and call set_time for every field in the field list.
     */
    bool set_time(double time);

    /**
     * If the field returns a FieldEnum and is constant on the given region, the method return true and
     * set @p value to the constant value on the given region. Otherwise (non constant field, other return type) it returns false.
     *
     * TODO: replace with more general method get_const_value
     */
    bool get_constant_enum_value(RegionIdx r_idx,  FieldEnum &value) const;

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
    virtual typename Value::return_type const &value(const Point<spacedim> &p, const ElementAccessor<spacedim> &elm);

    /**
     * Returns std::vector of scalar values in several points at once. The base class implements
     * trivial implementation using the @p value(,,) method. This is not optimal as it involves lot of virtual calls,
     * but this overhead can be negligible for more complex fields as Python of Formula.
     */
    virtual void value_list(const std::vector< Point<spacedim> >  &point_list, const  ElementAccessor<spacedim> &elm,
                       std::vector<typename Value::return_type>  &value_list);

private:
    /**
     *  Check that whole field list (@p region_fields_) is set, possibly use default values for unset regions.
     */
    void check_initialized_region_fields_();

    /**
     * If this pointer is set, turn off check of initialization in the set_time method on the regions
     * where the method get_constant_enum_value of the control field returns value from @p no_check_values_.
     */
    const Field<spacedim, typename FieldValue<spacedim>::Enum > *no_check_control_field_;
    std::vector<FieldEnum> no_check_values_;

    /**
     * Table with pointers to fields on individual regions.
     */
    std::vector< boost::shared_ptr< FieldBaseType > > region_fields_;

    /**
     * True after check_initialized_region_fields_ is called. That happen at first call of the set_time method.
     */
    bool is_fully_initialized_;
};



/**
 * Same as Field<...> but for boundary regions.
 */
template<int spacedim, class Value>
class BCField : public Field<spacedim, Value> {
public:
    BCField();
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
    virtual IT::AbstractRecord &get_input_type();

    virtual IT::AbstractRecord make_input_tree();

    /**
     * Abstract method for initialization of the field on one region.
     */
    virtual void set_from_input(const RegionSet &domain, const Input::AbstractRecord &rec);

    /**
     * Abstract method to update field to the new time level.
     * Implemented by in class template Field<...>.
     *
     * Return true if the value of the field was changed on some region.
     * The returned value is also stored in @p changed_during_set_time data member.
     */
    virtual bool set_time(double time);

    /**
     * We have to override the @p set_mesh method in order to call set_mesh method for subfields.
     */
    virtual void set_mesh(Mesh *mesh);

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
inline typename Value::return_type const & Field<spacedim,Value>::value(const Point<spacedim> &p, const ElementAccessor<spacedim> &elm)  {
    ASSERT(elm.region_idx().idx() < region_fields_.size(), "Region idx %u out of range %lu, field: %s\n",
           elm.region_idx().idx(), (unsigned long int) region_fields_.size(), this->name_.c_str());
    ASSERT( region_fields_[elm.region_idx().idx()] , "Null field ptr on region id: %d, field: %s\n", elm.region().id(), this->name_.c_str());
    return region_fields_[elm.region_idx().idx()]->value(p,elm);
}



template<int spacedim, class Value>
inline void Field<spacedim,Value>::value_list(const std::vector< Point<spacedim> >  &point_list, const ElementAccessor<spacedim> &elm,
                   std::vector<typename Value::return_type>  &value_list)
{
    region_fields_[elm.region_idx().idx()]->value_list(point_list,elm, value_list);
}





/****************************************************************************
 *  Macros for explicit instantiation of particular field class template.
 */


// Instantiation of fields with values dependent of the dimension of range space
#define INSTANCE_DIM_DEP_VALUES( field, dim_from, dim_to)                                                               \
template class field<dim_from, FieldValue<dim_to>::VectorFixed >;                       \
template class field<dim_from, FieldValue<dim_to>::TensorFixed >;                       \

// Instantiation of fields with domain in the ambient space of dimension @p dim_from
#define INSTANCE_TO_ALL(field, dim_from) \
template class field<dim_from, FieldValue<0>::Enum >;                       \
template class field<dim_from, FieldValue<0>::EnumVector >;                \
template class field<dim_from, FieldValue<0>::Integer >;                       \
template class field<dim_from, FieldValue<0>::Scalar >;                       \
template class field<dim_from, FieldValue<0>::Vector >;                         \
\
INSTANCE_DIM_DEP_VALUES( field, dim_from, 2) \
INSTANCE_DIM_DEP_VALUES( field, dim_from, 3) \

/*#define INSTANCE_ALL(field) \
INSTANCE_TO_ALL(field, 0) \
INSTANCE_TO_ALL( field, 1) \
INSTANCE_TO_ALL( field, 2) \
INSTANCE_TO_ALL( field, 3) */

// All instances of one field class template @p field.
// currently we need only fields on 3D ambient space (and 2D for some tests)
// so this is to save compilation time and avoid memory problems on the test server
#define INSTANCE_ALL(field) \
INSTANCE_TO_ALL( field, 2)  \
INSTANCE_TO_ALL( field, 3)

#endif /* FUNCTION_BASE_HH_ */
