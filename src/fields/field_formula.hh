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


#include "system/system.hh"
#include "fields/field_algo_base.hh"
#include "mesh/point.hh"
#include "input/factory.hh"

#include <string>
using namespace std;

class FunctionParser;

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

    virtual void init_from_input(const Input::Record &rec);

    /**
     * For time dependent formulas returns always true. For time independent formulas returns true only for the first time.
     */
    bool set_time(const TimeStep &time) override;

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
    // FieldValue_ wrapper for matrix of strings
    typedef FieldValue_<Value::NRows_, Value::NCols_, std::string> StringValue;

    // StringValue::return_type == StringTensor, which behaves like arma::mat<string>
    typename StringValue::return_type formula_matrix_;

    // FieldValue_ wrapper for unified reading of the input
    StringValue formula_matrix_helper_;

    // Matrix of parsers corresponding to the formula matrix returned by formula_matrix_helper_
    std::vector< std::vector<FunctionParser> > parser_matrix_;

    // Full address of the FiledFormula 'value' key.
    // Necessary in the case of an error during parsing.
    std::string value_input_address_;

    /// Registrar of class to factory
    static const int registrar;


};



#endif /* FIELD_FORMULA_HH_ */
