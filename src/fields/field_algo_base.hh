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
 * @file    field_algo_base.hh
 * @brief   
 * @todo
 * - better tests:
 *   - common set of quantities with different kind of values (scalar, vector, tensor, discrete, ..),
 *     common points and elements for evaluation
 *   - for individual Field implementations have:
 *     - different input
 *     - possibly different EPETCT_EQ tests, but rather have majority common
 */

#ifndef field_algo_base_HH_
#define field_algo_base_HH_

#include <string>
#include <memory>

#include <boost/type_traits.hpp>

#include "input/input_type_forward.hh"
#include "input/accessors_forward.hh"

#include "mesh/accessors.hh"
#include "mesh/point.hh"
#include "fields/field_values.hh"
#include "tools/time_governor.hh"



/**
 * Indication of special field states. Returned by Field<>::field_result.
 * Individual states have values corresponding to week ordering of the states according
 * to the exactness of the value. May possibly be helpful in implementation, e.g.
 * one can use (field_result >= result_constant) to check that the field is constant on given region.
 */
typedef enum  {
    result_none=0,      // field not set
    result_other=1,     // field initialized but no particular result information
    result_constant=2,  // spatially constant result
    result_zeros=10,    // zero scalar, vector, or tensor
    result_ones=20,     // all elements equal to 1.0
    result_eye=21       // identity tensor

} FieldResult;




/**
 * Base class for space-time function classes.
 */
template <int spacedim, class Value>
class FieldAlgorithmBase {
public:
       // expose template parameters
       typedef typename Space<spacedim>::Point Point;
       static const unsigned int spacedim_=spacedim;
       static constexpr bool is_enum_valued = boost::is_same<typename Value::element_type, FieldEnum>::value;


       /**
        * Kind of default constructor , with possible setting of the initial time.
        * Fields that returns variable size vectors accepts number of components @p n_comp.
        */
       FieldAlgorithmBase(unsigned int n_comp=0);

       /**
        * Returns template parameters as string in order to distinguish name of Abstracts
        * for initialization of different instances of the FieldBase template.
        */
       static std::string template_name();

       /**
        * Returns whole tree of input types for FieldBase with all descendants based on element input type (namely for FieldConstant)
        * given by element_input_type pointer.
        */
       static Input::Type::Abstract & get_input_type();

       /**
        * Returns parameterized whole tree of input types for FieldBase with all descendants based on element input type (namely
        * for FieldConstant) given by element_input_type pointer.
        */
       static const Input::Type::Instance & get_input_type_instance( Input::Type::Selection value_selection=Input::Type::Selection() );

       /**
        * This static method gets accessor to abstract record with function input,
        * dispatch to correct constructor and initialize appropriate function object from the input.
        * Returns shared pointer to  FunctionBase<>.
        */
       static std::shared_ptr< FieldAlgorithmBase<spacedim, Value> >
           function_factory(const Input::AbstractRecord &rec, unsigned int n_comp=0);

       /**
        *  Function can provide way to initialize itself from the input data.
        *
        *  TODO: make protected, should be called through function factory
        */
       virtual void init_from_input(const Input::Record &rec);

       /**
        * Set new time value. Some Fields may and some may not implement time dependent values and
        * possibly various types of interpolation. There can not be unified approach to interpolation (at least not on this abstraction level)
        * since some fields (FieldFormula, FieldPython) provides naturally time dependent functions other fields like (FieldConstant, ...), however,
        * can be equipped by various time interpolation schemes. In future, we obviously need time interpolation of higher order so that
        * we can use ODE integrators of higher order.
        *
        * The method returns true if the value of the field has changed in the new time step.
        */
       virtual bool set_time(const TimeStep &time);

       /**
        * Is used only by some Field implementations, but can be used to check validity of incoming ElementAccessor in value methods.
        *
        * Optional parameter @p boundary_domain can be used to specify, that the field will be evaluated only on the boundary part of the mesh.
        * TODO: make separate mesh for the boundary, then we can drop this parameter.
        */
       virtual void set_mesh(const Mesh *mesh, bool boundary_domain);

       /**
        * Sets @p component_idx_
        */
       void set_component_idx(unsigned int idx)
       { this->component_idx_ = idx; }

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
       { ASSERT(false).error("Not implemented yet."); return 0.0; }

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
       virtual typename Value::return_type const &value(const Point &p, const ElementAccessor<spacedim> &elm)=0;

       /**
        * Returns std::vector of scalar values in several points at once. The base class implements
        * trivial implementation using the @p value(,,) method. This is not optimal as it involves lot of virtual calls,
        * but this overhead can be negligible for more complex fields as Python of Formula.
        *
        * FieldAlgorithmBase provides a slow implementation using the value() method. Derived Field can implement its value_list method
        * as call of FieldAlgoritmBase<...>::value_list().
        */
       virtual void value_list(const std::vector< Point >  &point_list, const ElementAccessor<spacedim> &elm,
                          std::vector<typename Value::return_type>  &value_list)=0;

       /**
        * Virtual destructor.
        */
       virtual ~FieldAlgorithmBase() {}


protected:
       /// Actual time level; initial value is -infinity.
       TimeStep time_;
       /// Last value, prevents passing large values (vectors) by value.
       Value value_;
       typename Value::return_type r_value_;
       /// Indicator of particular values (zero, one) constant over space.
       FieldResult field_result_;
       /// Specify if the field is part of a MultiField and which component it is
       unsigned int component_idx_;
};


#endif /* FUNCTION_BASE_HH_ */
