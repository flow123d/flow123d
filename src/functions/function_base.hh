/*
 * function_base.hh
 *
 *  Created on: Aug 31, 2012
 *      Author: jb
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
       static FunctionBase<dim> *function_factory(Input::AbstractRecord &rec);

       /**
        *  Function can provide way to initialize itself from the input data.
        */
       virtual void init_from_input(Input::Record &rec);

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
