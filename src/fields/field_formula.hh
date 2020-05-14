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
    virtual void value_list (const Armor::array &point_list, const ElementAccessor<spacedim> &elm,
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



#endif /* FIELD_FORMULA_HH_ */
