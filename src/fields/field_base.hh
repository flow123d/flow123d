/*
 * field_base.hh
 *
 *  Created on: Aug 31, 2012
 *      Author: jb
 */

/**
 * TODO:
 * - FunctionBase (as well as all functions) will be templated by the type of the returned value @p Val
 *   and dimension
 * - methods:
 *   /// returns the value (for nontrivial Values this involves copy constructor)
 *   virtual Val value(Point<spacedim>, ElementAccessor<dim,spacedim>)
 *   /// Returns value through reference, the returned ResultType indicate zero, identity, not def and possibly other
 *   /// particular values. For complex 'Val' the values are not filled for nontrivial ResultType, i.e. we assume that
 *   /// there is an check of these particular cases. We may provide default resolution function.
 *   virtual ResultType value(Point<spacedim>, ElementAccessor<dim,spacedim>, Val &val);
 *   virtual void value_list(std::vector<Point<spacedim> >, ElementAccessor<dim,spacedim>, std::vector<Val> &, std::vecto<ResultType>& );
 *
 * - Question: how to treat parameter <dim> of ElementAccessors
 *   What we use from ElementAccessor?
 *   1) material number
 *   2) access to data on the same or refined mesh, i.e. make DoFAccessor from it
 *      identification of mesh, submesh, level, index in level
 *   3) Coordinates to interpolate from different mesh
 *
 *   Seems that nothing depends on <dim>
 *
 */

#include <fields/field_all.hh>

#ifndef FIELD_BASE_HH_
#define FIELD_BASE_HH_

#include <string>

#include "input/input_type.hh"
#include "input/accessors.hh"

#include "mesh/accessors.hh"
#include "mesh/point.hh"

#include "fields/field_values.hh"


/// Result type have sense only for larger Value types like vectors and tensors.
typedef enum  {
    result_zero,
    result_one,
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
       FieldBase(const double init_time=0.0, unsigned int n_comp=0);

       /**
        * Returns template parameters as string in order to distinguish Fields Input::Type name.
        */
       static std::string template_name();

       /**
        * Returns input type specification. As static function this depends on the template parameters.
        */
       static Input::Type::AbstractRecord &get_input_type();

       /**
        * This static method gets accessor to abstract record with function input,
        * dispatch to correct constructor and initialize appropriate function object from the input.
        * Returns pointer to  FunctionBase<>.
        */
       static FieldBase<spacedim, Value> *function_factory(Input::AbstractRecord rec,
               double init_time=0.0, unsigned int n_comp=0);

       /**
        *  Function can provide way to initialize itself from the input data.
        */
       virtual void init_from_input(Input::Record rec);

       /**
        * Set new time value.
        */
       virtual void set_time(double time);

       /**
        * Method for getting some information about next time where the function change its character.
        * Used to add appropriate TimeMarks.
        * TODO: think what kind of information we may need, is the next time value enough?
        */
       virtual double next_change_time()
       { ASSERT(0, "Not implemented yet."); }

       /**
        * Returns one value in one given point @p on an element given by ElementAccessor @p elm.
        * It returns reference to he actual value in order to avoid temporaries for vector and tensor values.
        *
        * This method just call the later one @p value(Point, ElementAccessor, Value) and drops the FieldResult.
        */
       virtual typename Value::return_type &value(const Point<spacedim> &p, ElementAccessor<spacedim> &elm);

       /**
        * Pure virtual method. At least this has to be implemented by descendants.
        * Returns one value in one given point. ResultType can be used to avoid some costly calculation if the result is trivial.
        */
       virtual FieldResult value(const Point<spacedim> &p, ElementAccessor<spacedim> &elm, typename Value::return_type &value) =0;

       /**
        * Returns std::vector of scalar values in several points at once. The base class implements
        * trivial implementation using the @p value(,,) method. This is not optimal as it involves lot of virtual calls,
        * but this overhead can be negligable for more complex fields as Python of Formula.
        */
       virtual void value_list(const std::vector< Point<spacedim> >  &point_list, ElementAccessor<spacedim> &elm,
                          std::vector<typename Value::return_type>  &value_list,
                          std::vector<FieldResult> &result_list);


protected:
       /// Actual time level
       double time_;
       /// Last value, prevents passing large values (vectors) by value.
       Value value_;
       typename Value::return_type r_value_;
};




#endif /* FUNCTION_BASE_HH_ */
