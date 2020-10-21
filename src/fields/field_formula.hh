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


#define BOOST_MATH_DISABLE_FLOAT128


#include <stdio.h>                      // for sprintf
#include <boost/exception/info.hpp>     // for operator<<, error_info::error...
#include <string>                       // for operator==, string
#include <vector>                       // for vector
#include <memory>
#include <armadillo>
#include "fields/field_algo_base.hh"    // for FieldAlgorithmBase
#include "fields/field_values.hh"       // for FieldValue<>::Enum, FieldValu...
#include "fields/field_set.hh"
#include "input/accessors.hh"           // for ExcAccessorForNullStorage
#include "input/accessors_impl.hh"      // for Record::val
#include "input/storage.hh"             // for ExcStorageTypeMismatch
#include "input/type_record.hh"         // for Record::ExcRecordKeyNotFound
#include "system/exceptions.hh"         // for ExcAssertMsg::~ExcAssertMsg
#include "tools/time_governor.hh"       // for TimeStep
#include "include/assert.hh"            // bparser
#include "include/parser.hh"            // bparser

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

    /// Size of data processed in BParser.
    static constexpr unsigned int bparser_vec_size = 128;

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

    void cache_update(FieldValueCache<typename Value::element_type> &data_cache,
			ElementCacheMap &cache_map, unsigned int region_idx) override;

    /**
     * Overload @p FieldAlgorithmBase::cache_reinit
     *
     * Reinit Bparser::ArenaAlloc data member.
     */
    void cache_reinit(const ElementCacheMap &cache_map) override;


    /**
     * Set reference of FieldSet.
     */
    std::vector<const FieldCommon *> set_dependency(FieldSet &field_set) override;

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
    std::vector< bparser::Parser > b_parser_;

    /// Accessor to Input::Record
    Input::Record in_rec_;

    /// Surface depth object calculate distance from surface.
    std::shared_ptr<SurfaceDepth> surface_depth_;

    /// Flag indicates if depth variable 'd' is used in formula
    bool has_depth_var_;

    /// Flag indicates first call of set_time method, when FunctionParsers in parser_matrix_ must be initialized
    bool first_time_set_;

    /// Holds FieldSet, allows evaluate values of Fields in formula expressions.
    FieldSet *field_set_;

    /// Arena object providing data arrays
    bparser::ArenaAlloc * arena_alloc_;

    // BParser data arrays
	double *X_;     ///< Data of coordinates, holds blocks of x, y, z
	double *x_;     ///< Coordinates x, part of previous array
	double *y_;     ///< Coordinates y, part of previous array
	double *z_;     ///< Coordinates z, part of previous array
	double *d_;     ///< Surface depth variable, used optionally if 'd' variable is set
	double *res_;   ///< Result vector of BParser
	uint *subsets_; ///< Subsets indices in range 0 ... n-1

    /// Registrar of class to factory
    static const int registrar;


};

// Necessary to linking.
template <int spacedim, class Value>
constexpr unsigned int FieldFormula<spacedim, Value>::bparser_vec_size;



#endif /* FIELD_FORMULA_HH_ */
