/*
 * function_base.hh
 *
 *  Created on: Aug 31, 2012
 *      Author: jb
 */

#ifndef FUNCTION_BASE_HH_
#define FUNCTION_BASE_HH_

#include <armadillo>

template <int dim>
class FunctionBase {
public:
       typedef arma::vec::fixed<dim> Point;

       FunctionBase(const unsigned int n_components=1, const double init_time=0.0)
       : n_components_(n_components), time_(init_time)
       {}

       static Input::Type::AbstractRecord &get_input_type();

       virtual void init_from_input(Input::Record &in_rec) = 0;

       /**
        * Set new time value.
        */
       virtual void set_time(const double time)
       { time_ = time; }

       /**
        * Returns one scalar value in one given point.
        */
       virtual double value(const Point &p, const unsigned int  component = 0) const =0;

       /**
        * Returns one vector value in one given point.
        */
       virtual void   vector_value(const Point &p, std::vector<double>     &value) const =0;

       /**
        * Returns std::vector of scalar values in several points at once.
        */
       virtual void   value_list (const std::vector< Point >  &point_list,
                          std::vector<double>         &value_list,
                          const unsigned int  component = 0) const =0;

       /**
        * Returns std::vector of scalar values in several points at once.
        */
       virtual void   vector_value_list (const std::vector< Point >    &point_list,
                                 std::vector< std::vector<double> >      &value_list) const=0;

       virtual ~FunctionBase() {}

protected:
       unsigned int n_components_;
       double time_;
};



#endif /* FUNCTION_BASE_HH_ */
