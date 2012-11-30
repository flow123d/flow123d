/*
 * function_base.hh
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


#include <functions/functions_all.hh>

#ifndef FUNCTION_BASE_HH_
#define FUNCTION_BASE_HH_

#include <armadillo>
#include "input/input_type.hh"
#include "mesh/elements.h"



template<int dim>
class FunctionBase;

/**
 * Base class for space-time function classes.
 */
template <int dim>
class FunctionBase {
public:
       typedef arma::vec::fixed<dim> Point;

       FunctionBase(const unsigned int n_components=1, const double init_time=0.0);

       static Input::Type::AbstractRecord &get_input_type();

       /**
        * This static method gets accessor to abstract record with function input,
        * dispatch to correct constructor and initialize appropriate function object from the input.
        * Returns pointer to general FunctionBase.
        */
       static FunctionBase<dim> *function_factory(Input::AbstractRecord rec, const unsigned int n_comp=0, const double init_time=0.0);

       /**
        *  Function can provide way to initialize itself from the input data.
        */
       virtual void init_from_input(Input::Record rec);

       /**
        * Set new time value.
        */
       virtual void set_time(double time);

       /**
        * Set element for interpolation
        */
       virtual void set_element(Element *element);


       /**
        * Returns one scalar value in one given point.
        */
       virtual double value(const Point &p, const unsigned int  component = 0) const =0;

       /**
        * Returns one vector value in one given point.
        */
       virtual void vector_value(const Point &p, std::vector<double>     &value) const =0;

       /**
        * Returns std::vector of scalar values in several points at once.
        */
       virtual void value_list (const std::vector< Point >  &point_list,
                          std::vector<double>         &value_list,
                          const unsigned int  component = 0) const =0;

       /**
        * Returns std::vector of scalar values in several points at once.
        */
       virtual void vector_value_list (const std::vector< Point >    &point_list,
                                 std::vector< std::vector<double> >      &value_list) const=0;


protected:
       unsigned int n_components_;
       double time_;
       Element *element_;

};




#endif /* FUNCTION_BASE_HH_ */
