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
 * @file    field_formula.hh
 * @brief   
 */

#ifndef FIELD_FORMULA_HH_
#define FIELD_FORMULA_HH_


#include <stdio.h>                      // for sprintf
#include <boost/exception/info.hpp>     // for operator<<, error_info::error...
#include <string>                       // for operator==, string
#include <vector>                       // for vector
#include <memory>
#include <armadillo>
#include "fields/field_algo_base.hh"    // for FieldAlgorithmBase
#include "fields/field_values.hh"       // for FieldValue<>::Enum, FieldValu...
#include "input/accessors.hh"           // for ExcAccessorForNullStorage
#include "input/accessors_impl.hh"      // for Record::val
#include "input/storage.hh"             // for ExcStorageTypeMismatch
#include "input/type_record.hh"         // for Record::ExcRecordKeyNotFound
#include "system/exceptions.hh"         // for ExcAssertMsg::~ExcAssertMsg
#include "tools/time_governor.hh"       // for TimeStep

class FunctionParser;
template <int spacedim> class ElementAccessor;
class SurfaceDepth;

using namespace std;

/**
 * Class representing fields given by runtime parsed formulas.
 *
 * Using library:
 * http://warp.povusers.org/FunctionParser/
 *
 * TODO:
 * correct support for discrete functions (use integer parser), actually we just convert double to int
 *
 */
template <int spacedim, class Value>
class FieldFormula : public FieldAlgorithmBase<spacedim, Value>
{
public:
    typedef typename FieldAlgorithmBase<spacedim, Value>::Point Point;
    typedef FieldAlgorithmBase<spacedim, Value> FactoryBaseType;

    FieldFormula(unsigned int n_comp=0);


    static const Input::Type::Record & get_input_type();

    virtual void init_from_input(const Input::Record &rec, const struct FieldAlgoBaseInitData& init_data);

    /**
     * For time dependent formulas returns always true. For time independent formulas returns true only for the first time.
     */
    bool set_time(const TimeStep &time) override;

    /**
     * Create SurfaceDepth object if surface region is set.
     *
     * See also description of the FieldBase<...>::set_mesh.
     */
    void set_mesh(const Mesh *mesh, bool boundary_domain) override;

    /**
     * Returns one value in one given point. ResultType can be used to avoid some costly calculation if the result is trivial.
     */
    virtual typename Value::return_type const &value(const Point &p, const ElementAccessor<spacedim> &elm);

    /**
     * Returns std::vector of scalar values in several points at once.
     */
    virtual void value_list (const std::vector< Point >  &point_list, const ElementAccessor<spacedim> &elm,
                       std::vector<typename Value::return_type>  &value_list);


    virtual ~FieldFormula();

private:
    typedef StringTensorInput<Value::NRows_,Value::NCols_> STI;

    /**
     * Evaluate depth variable if it is contained in formula.
     *
     * Return arma vec of point coordinates extended by depth value (or zero if depth is not contained.
     */
    inline arma::vec eval_depth_var(const Point &p);

    // StringValue::return_type == StringTensor, which behaves like arma::mat<string>
    StringTensor formula_matrix_;

    // Matrix of parsers corresponding to the formula matrix returned by formula_matrix_helper_
    std::vector< std::vector<FunctionParser> > parser_matrix_;

    /// Accessor to Input::Record
    Input::Record in_rec_;

    /// Surface depth object calculate distance from surface.
    std::shared_ptr<SurfaceDepth> surface_depth_;

    /// Flag indicates if depth variable 'd' is used in formula
    bool has_depth_var_;

    /// Flag indicates first call of set_time method, when FunctionParsers in parser_matrix_ must be initialized
    bool first_time_set_;

    /// Registrar of class to factory
    static const int registrar;


};

/**
 * Documentation of exprtk parser:
 *
 * FieldFormula uses exprtk library for parsing expressions. For correct functionality,
 * you can use follow operators and functions of exprtk module:
 *
 * +-------------+-------------------------------------------------+-------------+
 * | OPERATOR    | DESCRIPTION,                                    | EXAMPLE     |
 * | OR FUNCTION | NOTE                                            | OF USAGE    |
 * +-------------+-------------------------------------------------+-------------+
 * | +           | Addition operator.                              | x + y       |
 * +-------------+-------------------------------------------------+-------------+
 * | -           | Subtraction operator.                           | x - y       |
 * +-------------+-------------------------------------------------+-------------+
 * | *           | Multiplication operator.                        | x * y       |
 * +-------------+-------------------------------------------------+-------------+
 * | /           | Division operator.                              | x / y       |
 * +-------------+-------------------------------------------------+-------------+
 * | %           | Modulus operator.                               | x % y       |
 * +-------------+-------------------------------------------------+-------------+
 * | ^           | Exponent operator (x to the power of y).        | x ^ y       |
 * +-------------+-------------------------------------------------+-------------+
 * | :=          | Assignment operator.                            | i := 0      |
 * +-------------+-------------------------------------------------+-------------+
 * | var         | Definition of variable.                         | var i       |
 * +-------------+-------------------------------------------------+-------------+
 * | ;           | Statement terminator, same meaning as in C++.   | var i := 5; |
 * |             | Note: Not obligatory for single statement.      | x + i       |
 * +-------------+-------------------------------------------------+-------------+
 * | == or =     | Comparison operator equal to.                   | x == y      |
 * +-------------+-------------------------------------------------+-------------+
 * | <> or !=    | Comparison operator not equal to.               | x <> y      |
 * |             |                                                 | x != y      |
 * +-------------+-------------------------------------------------+-------------+
 * | <           | Comparison operator less than.                  | x < y       |
 * +-------------+-------------------------------------------------+-------------+
 * | <=          | Comparison operator less than or equal.         | x <= y      |
 * +-------------+-------------------------------------------------+-------------+
 * | >           | Comparison operator greater than.               | x > y       |
 * +-------------+-------------------------------------------------+-------------+
 * | >=          | Comparison operator greater than or equal.      | x >= y      |
 * +-------------+-------------------------------------------------+-------------+
 * | abs         | Absolute value function.                        | abs(x)      |
 * +-------------+-------------------------------------------------+-------------+
 * | exp         | Exponentioal function (e to the power of x).    | exp(x)      |
 * +-------------+-------------------------------------------------+-------------+
 * | log         | Natural logarithm function.                     | log(x)      |
 * +-------------+-------------------------------------------------+-------------+
 * | log10       | Base 10 logarithm function.                     | log10(x)    |
 * +-------------+-------------------------------------------------+-------------+
 * | log2        | Base 2 logarithm function.                      | log2(x)     |
 * +-------------+-------------------------------------------------+-------------+
 * | logN        | Base N logarithm function.                      | logN(x,8)   |
 * |             | Note: N is positive integer.                    |             |
 * +-------------+-------------------------------------------------+-------------+
 * | sin         | Sine function.                                  | sin(x)      |
 * |             | Note: Equivalent syntax of cos, tan, cot.       |             |
 * +-------------+-------------------------------------------------+-------------+
 * | asin        | Arc sine function of expression in radians.     | asin(x)     |
 * |             | Note: Equivalent syntax of acos, atan.          |             |
 * +-------------+-------------------------------------------------+-------------+
 * | min         | Smallest value of all the inputs.               | min(x,y,z)  |
 * +-------------+-------------------------------------------------+-------------+
 * | max         | Largest value of all the inputs.                | max(x,y,z)  |
 * +-------------+-------------------------------------------------+-------------+
 * | ? :         | Ternary conditional statement, similar to C++   | x>0 ? y : z |
 * |             | ternary operator.                               |             |
 * +-------------+-------------------------------------------------+-------------+
 *
 *
 * Speed of evaluation of selected functions and operators:
 * If we assume that the evaluation of the addition operation takes 1 TUNIT, the following
 * table shows how many times is (approximately) the evaluation of each function slower.
 *
 * +------------------------+-------+
 * | FUNCTIONS              | TUNIT |
 * +------------------------+-------+
 * | min, max               |     1 |
 * +------------------------+-------+
 * | abs                    |     1 |
 * +------------------------+-------+
 * | ? :                    |     2 |
 * +------------------------+-------+
 * | trigonometry functions |    50 |
 * +------------------------+-------+
 * | logarithmic functions  |    55 |
 * +------------------------+-------+
 * | exp                    |    60 |
 * +------------------------+-------+
 * | ^                      |   130 |
 * +------------------------+-------+
 *
 *
 * Notes for developers:
 * The parser does not throw an exception if expression is incorrect. We must check
 * return value of method exprtk::parser<Type>::compile() and if it returns 'false'
 * we must throw a suitable exception.
 *
 * Full example of parser code with using vector view:
@code
typedef exprtk::symbol_table<double> symbol_table_t;
typedef exprtk::expression<double>     expression_t;
typedef exprtk::parser<double>             parser_t;

// Define input vectors of x, y, z data and output result vector (all vectors with same sizes)
std::vector<double> x_vec = ...
std::vector<double> y_vec = ...
std::vector<double> z_vec = ...
std::vector<double> result(0.0, x_vec.size());

// Evaluated expression
std::string expression_string = " result_vec := x + y + z ";

// Definitions of vector views
exprtk::vector_view<double> x_view = exprtk::make_vector_view(x_vec,x_vec.size());
exprtk::vector_view<double> y_view = exprtk::make_vector_view(y_vec,y_vec.size());
exprtk::vector_view<double> z_view = exprtk::make_vector_view(z_vec,z_vec.size());
exprtk::vector_view<double> r_view = exprtk::make_vector_view(result,result.size());

// Definition of symbol table, needs method add constant for declaration (pi, e ...)
symbol_table_t symbol_table;
symbol_table.add_vector("x",x_view);
symbol_table.add_vector("y",y_view);
symbol_table.add_vector("z",z_view);
symbol_table.add_vector("result_vec",r_view);
symbol_table.add_constants();

// Compilation and parsing
expression_t expression;
expression.register_symbol_table(symbol_table);
parser_t parser;
parser.compile(expression_string,expression); // Needs check correct expression
expression.value();

// Processing of result
for (double val : result)
    some_process(val);
@endcode
 *
We can use method 'rebase' of all vector views and call 'value' repeatedly with
different input vectors. Different parts of code are:
@code
// Define set of input vectors and result vector with same sies
std::vector< std::vector<double> > x_vecs = ...
std::vector< std::vector<double> > y_vecs = ...
std::vector< std::vector<double> > z_vecs = ...
std::vector<double> result(0.0, x_vecs[0].size());

...
// Parsing and processing of result
for (unsigned int i=0; i<x_vecs.size(); ++i)
{
    x_view.rebase(x_vecs[i].data()); // update vectors
    z_view.rebase(y_vecs[i].data());
    y_view.rebase(z_vecs[i].data());
    expression.value();
    // processing of result
}
@endcode
 *
 */

#endif /* FIELD_FORMULA_HH_ */
