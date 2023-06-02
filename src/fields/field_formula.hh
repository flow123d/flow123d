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
#include <string>                       // for operator==, string
#include <vector>                       // for vector
#include <memory>
#include <armadillo>
#include <map>
#include "fields/field_algo_base.hh"    // for FieldAlgorithmBase
#include "fields/field_values.hh"       // for FieldValue<>::Enum, FieldValu...
#include "fields/field_set.hh"
#include "input/accessors.hh"           // for ExcAccessorForNullStorage
#include "input/accessors_impl.hh"      // for Record::val
#include "input/storage.hh"             // for ExcStorageTypeMismatch
#include "input/type_record.hh"         // for Record::ExcRecordKeyNotFound
#include "input/input_exception.hh"     // for ExcAssertMsg::~ExcAssertMsg
#include "system/exceptions.hh"         // for ExcAssertMsg::~ExcAssertMsg
#include "tools/time_governor.hh"       // for TimeStep
//#include "include/assert.hh"            // bparser
#include "parser.hh"            // bparser

class FunctionParser;
template <int spacedim> class ElementAccessor;
class SurfaceDepth;

using namespace std;

/**
 * Class representing fields given by runtime parsed formulas.
 *
 * Using libraries:
 * https://github.com/flow123d/bparser/
 * http://warp.povusers.org/FunctionParser/ (gradually replaced by BParser)
 *
 * Allows parsing:
 * - base variables: coordinates (x,y,z), time (t), surface depth (d); constants: e, pi
 * - standard functions: trigonometric functions, min and max function, exponential function, ternary operator etc.
 * - expressions dependendent on other fields
 */
template <int spacedim, class Value>
class FieldFormula : public FieldAlgorithmBase<spacedim, Value>
{
public:
    typedef typename FieldAlgorithmBase<spacedim, Value>::Point Point;
    typedef FieldAlgorithmBase<spacedim, Value> FactoryBaseType;

    TYPEDEF_ERR_INFO(EI_Field, std::string);
    DECLARE_INPUT_EXCEPTION(ExcNotDoubleField,
            << "Can not use integer valued field " << EI_Field::qval << " in the formula: \n");

    TYPEDEF_ERR_INFO(EI_BParserMsg, std::string);
    TYPEDEF_ERR_INFO(EI_Formula, std::string);
    DECLARE_INPUT_EXCEPTION(ExcParserError,
            << "Parsing in " << EI_BParserMsg::val << " in the formula: " << EI_Formula::qval << "\n");

    // Temporary exception of FParser. TODO remove at the same time as FParser
    TYPEDEF_ERR_INFO(EI_FParserMsg, std::string);
    TYPEDEF_ERR_INFO(EI_Row, unsigned int);
    TYPEDEF_ERR_INFO(EI_Col, unsigned int);
    DECLARE_INPUT_EXCEPTION(ExcFParserError,
            << "ParserError: " << EI_FParserMsg::val << "\n in the FieldFormula[" << EI_Row::val
			<< "][" << EI_Row::val << "] == " << EI_Formula::qval << " \n");

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

    void cache_update(FieldValueCache<typename Value::element_type> &data_cache,
			ElementCacheMap &cache_map, unsigned int region_patch_idx) override;

    /**
     * Set reference of FieldSet.
     */
    std::vector<const FieldCommon *> set_dependency(FieldSet &field_set) override;

    /**
     * Overload @p FieldAlgorithmBase::cache_reinit
     *
     * Reinit arena data member.
     */
    void cache_reinit(const ElementCacheMap &cache_map) override;

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

    /// Flag indicates if depth variable 'd' is used in formula - obsolete parameter of FParser
    bool has_depth_var_;

    /// Flag indicates if time variable 't' is used in formula - parameter of BParser
    bool has_time_;

    /// Helper variable for construct of arena, holds sum of sizes (over shape) of all dependent fields.
    uint sum_shape_sizes_;

    /// Flag indicates first call of set_time method, when FunctionParsers in parser_matrix_ must be initialized
    bool first_time_set_;

    /// Arena object providing data arrays
    bparser::ArenaAlloc * arena_alloc_;

    // BParser data arrays and variables
	double *x_;       ///< Coordinates x, part of previous array
	double *y_;       ///< Coordinates y, part of previous array
	double *z_;       ///< Coordinates z, part of previous array
	double *d_;       ///< Surface depth variable, used optionally if 'd' variable is set
	double *res_;     ///< Result vector of BParser
	uint *subsets_;   ///< Subsets indices in range 0 ... n-1
	std::vector<const FieldCommon * > required_fields_;

	/**
	 * Data of fields evaluated in expressions.
	 *
	 * Temporary data member, we need to copy data from FieldValueCaches to arrays allocated in arena.
	 */
	std::unordered_map<const FieldCommon *, double *> eval_field_data_;

    /// Registrar of class to factory
    static const int registrar;


};


#endif /* FIELD_FORMULA_HH_ */
