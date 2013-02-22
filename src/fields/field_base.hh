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
        */
       virtual void set_time(double time);

       /**
        * Is used only by some Field imlementations, but can be used to check validity of incomming ElementAccessor in value methods.
        */
       virtual void set_mesh(Mesh *mesh);

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
       /// Actual time level
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
     * Setters. We need to store these information into the Field itself during construction of an EqData class in order to
     * allow creation of Input::Type tree and initialization of all fields in the EqData class by generic functions.
     */
    /// Set name of the field, used for naming particular key in EqData record.
    void set_name(const string & name);
    /// Set description of the field, used for description of corresponding key.
    void set_desc(const string & desc);
    void set_default(const IT::Default &dflt);

    /// Set number of components for run-time sized vectors.
    void set_n_comp( unsigned int n_comp);
    /// For Fields returning "Enum", we have to pass in corresponding Selection object.
    void set_selection( Input::Type::Selection *element_selection);
    void set_mesh(Mesh *mesh);

    /**
     * Getters.
     */
    const std::string &name() const;
    const std::string &desc() const;
    const IT::Default &get_default() const;
    bool is_bc() const;
    bool is_enum_valued() const;
    unsigned int n_comp() const;
    Mesh * mesh() const;

    /**
     * Returns input type of particular field instance, this is usually static member input_type of the corresponding FieldBase class (
     * with same template parameters), however, for fields returning "Enum" we have to create whole unique Input::Type hierarchy for
     * every instance since every such field use different Selection for initialization, even if all returns just unsigned int.
     */
    virtual IT::AbstractRecord &get_input_type() =0;

    virtual IT::AbstractRecord make_input_tree() =0;

    /**
     * Abstract method for initialization of the field on one region.
     */
    virtual void set_from_input(const RegionSet &domain, const Input::AbstractRecord &rec) =0;

    /**
     * Abstract method to update field to the new time.
     */
    virtual void set_time(double time) =0;

    /**
     * Virtual destructor.
     */
    virtual ~FieldCommonBase();

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
};





///Helper function.
template <class FieldBaseType>
IT::AbstractRecord get_input_type_resolution(
        Input::Type::Selection *sel,  const boost::true_type&)
{
    ASSERT( sel, "NULL pointer to selection in Field::get_input_type(), while Value==FieldEnum.\n");
    // create nonlocal copy that live as long as the object instance
    //rec = boost::make_shared<IT::AbstractRecord>(FieldBaseType::get_input_type(sel));
    return FieldBaseType::get_input_type(sel);
}
template <class FieldBaseType>
IT::AbstractRecord get_input_type_resolution(
        Input::Type::Selection *sel,  const boost::false_type&)
{
    return FieldBaseType::get_input_type(NULL);
}





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
 * TODO: simplify initialization of the class, currently it needs:
 * - default constructor .. OK
 * - set_mesh, possibly other FieldCommonBase setters
 * - optional set_from_input or set_field
 * - set_time
 *
 * there should be clearly three distinguish states: after default construction, having mesh (with closed regiondb - allocate table), checked_table
 * in first set_time
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
     * Returns input type of particular field instance, this is usually static member input_type of the corresponding FieldBase class (
     * with same template parameters), however, for fields returning "Enum" we have to create whole unique Input::Type hierarchy for
     * every instance since every such field use different Selection for initialization, even if all returns just unsigned int.
     */
    IT::AbstractRecord &get_input_type() {
        return FieldBaseType::input_type;
    }

    IT::AbstractRecord make_input_tree() {
        return get_input_type_resolution<FieldBaseType>( this->element_selection_ ,boost::is_same<typename Value::element_type, FieldEnum>());
    }

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
    void set_time(double time);

    /**
     * If the field returns a FieldEnum and is constant on the given region, the method return true and
     * set @p value to the constant value on the given region. Otherwise (non constant field, other return type) it returns false.
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
     * If this pointer is set, turn off check of initialization in the set_time method on the regions
     * where the method get_constant_enum_value of the constrol field returns value from @p no_check_values_.
     */
    const Field<spacedim, typename FieldValue<spacedim>::Enum > *no_check_control_field_;
    std::vector<FieldEnum> no_check_values_;

    std::vector< boost::shared_ptr< FieldBaseType > > region_fields_;


};



/**
 * Same as Field<...> but for boundary regions.
 */
template<int spacedim, class Value>
class BCField : public Field<spacedim, Value> {
public:
    BCField() { this->bc_=true; }
};



#endif /* FUNCTION_BASE_HH_ */
