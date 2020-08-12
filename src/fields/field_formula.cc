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
 * @file    field_formula.cc
 * @brief   
 */



#include "fields/field_formula.hh"
#include "fields/field_instances.hh"	// for instantiation macros
#include "fields/surface_depth.hh"
#include "fparser.hh"
#include "input/input_type.hh"
#include "include/arena_alloc.hh"       // bparser
#include <boost/foreach.hpp>
#include <boost/algorithm/string/replace.hpp>


/// Implementation.

namespace it = Input::Type;

FLOW123D_FORCE_LINK_IN_CHILD(field_formula)


template <int spacedim, class Value>
const Input::Type::Record & FieldFormula<spacedim, Value>::get_input_type()
{

    return it::Record("FieldFormula", FieldAlgorithmBase<spacedim,Value>::template_name()+" Field given by runtime interpreted formula.")
            .derive_from(FieldAlgorithmBase<spacedim, Value>::get_input_type())
            .copy_keys(FieldAlgorithmBase<spacedim, Value>::get_field_algo_common_keys())
            .declare_key("value", STI::get_input_type() , it::Default::obligatory(),
                                        "String, array of strings, or matrix of strings with formulas for individual "
                                        "entries of scalar, vector, or tensor value respectively.\n"
                                        "For vector values, you can use just one string to enter homogeneous vector.\n"
                                        "For square (($N\\times N$))-matrix values, you can use:\n\n"
                                        " - array of strings of size (($N$)) to enter diagonal matrix\n"
                                        " - array of strings of size (($\\frac12N(N+1)$)) to enter symmetric matrix (upper triangle, row by row)\n"
                                        " - just one string to enter (spatially variable) multiple of the unit matrix.\n"
                                        "Formula can contain variables ```x,y,z,t,d``` and usual operators and functions." )
			//.declare_key("unit", FieldAlgorithmBase<spacedim, Value>::get_input_type_unit_si(), it::Default::optional(),
			//							"Definition of unit.")
			.declare_key("surface_direction", it::String(), it::Default("\"0 0 1\""),
										"The vector used to project evaluation point onto the surface.")
			.declare_key("surface_region", it::String(), it::Default::optional(),
										"The name of region set considered as the surface. You have to set surface region if you "
										"want to use formula variable ```d```.")
	        .allow_auto_conversion("value")
			.close();
}



template <int spacedim, class Value>
const int FieldFormula<spacedim, Value>::registrar =
		Input::register_class< FieldFormula<spacedim, Value>, unsigned int >("FieldFormula") +
		FieldFormula<spacedim, Value>::get_input_type().size();



template <int spacedim, class Value>
FieldFormula<spacedim, Value>::FieldFormula( unsigned int n_comp)
: FieldAlgorithmBase<spacedim, Value>(n_comp),
  formula_matrix_(this->value_.n_rows(), this->value_.n_cols()),
  first_time_set_(true), field_set_(nullptr), arena_alloc_(nullptr)
{
	this->is_constant_in_space_ = false;
    parser_matrix_.resize(this->value_.n_rows());
    for(unsigned int row=0; row < this->value_.n_rows(); row++) {
        parser_matrix_[row].resize(this->value_.n_cols());
    }
    b_parser_.reserve(this->value_.n_rows()*this->value_.n_cols());
    for(unsigned int i=0; i < this->value_.n_rows()*this->value_.n_cols(); i++) {
        b_parser_.emplace_back( FieldFormula<spacedim, Value>::bparser_vec_size );
    }
}



template <int spacedim, class Value>
void FieldFormula<spacedim, Value>::init_from_input(const Input::Record &rec, const struct FieldAlgoBaseInitData& init_data) {
	this->init_unit_conversion_coefficient(rec, init_data);

	// read formulas form input
    STI::init_from_input( formula_matrix_, rec.val<typename STI::AccessType>("value") );
    in_rec_ = rec;
}


template <int spacedim, class Value>
bool FieldFormula<spacedim, Value>::set_time(const TimeStep &time) {


    bool any_parser_changed = false;
    std::string value_input_address = in_rec_.address_string();
    has_depth_var_ = false;
    this->is_constant_in_space_ = true; // set flag to true, then if found 'x', 'y', 'z' or 'd' reset to false


    std::string vars = string("x,y,z").substr(0, 2*spacedim-1);
    std::vector<bool> time_dependent(this->value_.n_rows() * this->value_.n_cols(), false);
    // update parsers
    for(unsigned int row=0; row < this->value_.n_rows(); row++)
        for(unsigned int col=0; col < this->value_.n_cols(); col++) {
            // get all variable names from the formula
            std::vector<std::string> var_list;

            FunctionParser tmp_parser;
            tmp_parser.AddConstant("Pi", 3.14159265358979323846);
            tmp_parser.AddConstant("E", 2.71828182845904523536);


#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wunused-variable"
            {
                int err=tmp_parser.ParseAndDeduceVariables(formula_matrix_.at(row,col), var_list);
                ASSERT(err != FunctionParser::FP_NO_ERROR)(tmp_parser.ErrorMsg()).error("ParseAndDeduceVariables error\n");
            }
#pragma GCC diagnostic pop

            BOOST_FOREACH(std::string &var_name, var_list ) {
                if (var_name == std::string("t") ) time_dependent[row*this->value_.n_rows()+col]=true;
                else if (var_name == std::string("d") ) {
                	this->is_constant_in_space_ = false;
                	if (surface_depth_)
                		has_depth_var_=true;
                	else
                    	WarningOut().fmt("Unset surface region. Variable '{}' in the FieldFormula[{}][{}] == '{}' will be set to zero\n at the input address:\n {} \n",
                                var_name, row, col, formula_matrix_.at(row,col), value_input_address );
                }
                else if (var_name == "x" || var_name == "y" || var_name == "z") {
                	this->is_constant_in_space_ = false;
                	continue;
                }
                else
                	WarningOut().fmt("Unknown variable '{}' in the  FieldFormula[{}][{}] == '{}'\n at the input address:\n {} \n",
                            var_name, row, col, formula_matrix_.at(row,col), value_input_address );
            }

            // Seems that we can not just add 't' constant to tmp_parser, since it was already Parsed.
            parser_matrix_[row][col].AddConstant("Pi", 3.14159265358979323846);
            parser_matrix_[row][col].AddConstant("E", 2.71828182845904523536);
            if (time_dependent[row*this->value_.n_rows()+col]) {
                parser_matrix_[row][col].AddConstant("t", time.end());
            }
        }

    if (has_depth_var_)
        vars += string(",d");

	// update parsers
	for(unsigned int row=0; row < this->value_.n_rows(); row++)
		for(unsigned int col=0; col < this->value_.n_cols(); col++) {
            // TODO:
            // - possibly add user defined constants and units here ...
            // - optimization; possibly parse only if time_dependent  || formula_matrix[][] has changed ...
            //parser_matrix_[row][col] = tmp_parser;
            if (time_dependent[row*this->value_.n_rows()+col] || first_time_set_ ) {
                parser_matrix_[row][col].Parse(formula_matrix_.at(row,col), vars);

                if ( parser_matrix_[row][col].GetParseErrorType() != FunctionParser::FP_NO_ERROR ) {
                    xprintf(UsrErr, "ParserError: %s\n in the FieldFormula[%d][%d] == '%s'\n at the input address:\n %s \n",
                        parser_matrix_[row][col].ErrorMsg(),
                        row,col,formula_matrix_.at(row,col).c_str(),
                        value_input_address.c_str());
                }

                parser_matrix_[row][col].Optimize();
                any_parser_changed = true;
            }


        }

    first_time_set_ = false;
    this->time_=time;
    return any_parser_changed;
}


template <int spacedim, class Value>
void FieldFormula<spacedim, Value>::set_mesh(const Mesh *mesh, FMT_UNUSED bool boundary_domain) {
    // create SurfaceDepth object if surface region is set
    std::string surface_region;
    if ( in_rec_.opt_val("surface_region", surface_region) ) {
        surface_depth_ = std::make_shared<SurfaceDepth>(mesh, surface_region, in_rec_.val<std::string>("surface_direction"));
    }
}


/**
 * Returns one value in one given point. ResultType can be used to avoid some costly calculation if the result is trivial.
 */
template <int spacedim, class Value>
typename Value::return_type const & FieldFormula<spacedim, Value>::value(const Point &p, FMT_UNUSED  const ElementAccessor<spacedim> &elm)
{

    auto p_depth = this->eval_depth_var(p);
    for(unsigned int row=0; row < this->value_.n_rows(); row++)
        for(unsigned int col=0; col < this->value_.n_cols(); col++) {
            this->value_(row,col) = this->unit_conversion_coefficient_ * parser_matrix_[row][col].Eval(p_depth.memptr());
        }
    return this->r_value_;
}


/**
 * Returns std::vector of scalar values in several points at once.
 */
template <int spacedim, class Value>
void FieldFormula<spacedim, Value>::value_list (const Armor::array &point_list, FMT_UNUSED const ElementAccessor<spacedim> &elm,
                   std::vector<typename Value::return_type>  &value_list)
{
	ASSERT_EQ( point_list.size(), value_list.size() );
    ASSERT_DBG( point_list.n_rows() == spacedim && point_list.n_cols() == 1).error("Invalid point size.\n");
    for(unsigned int i=0; i< point_list.size(); i++) {
        Value envelope(value_list[i]);
        ASSERT_EQ( envelope.n_rows(), this->value_.n_rows() )(i)(envelope.n_rows())(this->value_.n_rows())
        		.error("value_list['i'] has wrong number of rows\n");
        auto p_depth = this->eval_depth_var(point_list.vec<spacedim>(i));

        for(unsigned int row=0; row < this->value_.n_rows(); row++)
            for(unsigned int col=0; col < this->value_.n_cols(); col++) {
                envelope(row,col) = this->unit_conversion_coefficient_ * parser_matrix_[row][col].Eval(p_depth.memptr());
            }
    }
}


template <int spacedim, class Value>
void FieldFormula<spacedim, Value>::cache_update(FieldValueCache<typename Value::element_type> &data_cache,
        ElementCacheMap &cache_map, unsigned int region_idx)
{
    unsigned int reg_chunk_begin = cache_map.region_chunk_begin(region_idx);
    unsigned int reg_chunk_end = cache_map.region_chunk_end(region_idx);

    for (unsigned int i=reg_chunk_begin; i<reg_chunk_end; ++i) {
        // fill data vectors
    	x_[i] = field_set_->x().template mat<1,1>(i)(0);
        y_[i] = field_set_->y().template mat<1,1>(i)(0);
        z_[i] = field_set_->z().template mat<1,1>(i)(0);
        if (surface_depth_ && has_depth_var_) {
            Point p;
            p(0) = x_[i]; p(1) = y_[i]; p(2) = z_[i];
            // TODO computation of depth var needs better solution, probably we add field to FieldSet
            d_[i] = surface_depth_->compute_distance(p);
        }
        res_[i] = 0.0;
    }

    // Get vector of subsets as subarray
    uint subsets_begin = reg_chunk_begin / ElementCacheMap::simd_size_double;
    uint subsets_end = reg_chunk_end / ElementCacheMap::simd_size_double;
    std::vector<uint> subset_vec;
    subset_vec.assign(subsets_ + subsets_begin, subsets_ + subsets_end);

    for(unsigned int row=0; row < this->value_.n_rows(); row++)
        for(unsigned int col=0; col < this->value_.n_cols(); col++) {
            b_parser_[row*this->value_.n_cols()+col].set_subset(subset_vec);
            b_parser_[row*this->value_.n_cols()+col].run();
            for (unsigned int i=reg_chunk_begin; i<reg_chunk_end; ++i) {
                auto cache_val = data_cache.template mat<Value::NRows_, Value::NCols_>(i);
                cache_val(row, col) = res_[i];
                data_cache.set(i) = cache_val;
            }
        }
}


template <int spacedim, class Value>
void FieldFormula<spacedim, Value>::cache_reinit(const ElementCacheMap &cache_map)
{
	bool use_depth_var = (surface_depth_ && has_depth_var_); // TODO need better check of using 'd' variable (from 'Parser::variables()').
	if (arena_alloc_!=nullptr) {
	    delete arena_alloc_;
	}
	uint vec_size = cache_map.eval_points()->max_size() * ElementCacheMap::n_cached_elements;
	while (vec_size%ElementCacheMap::simd_size_double > 0) vec_size++; // alignment of block size
	// number of subset alignment to block size
	uint n_subsets = (vec_size+ElementCacheMap::simd_size_double-1) / ElementCacheMap::simd_size_double;
	uint n_vectors = (use_depth_var ? 5 : 4); // needs vectors of coordinates x, y, z, result vector and optionally d (depth)
	arena_alloc_ = new bparser::ArenaAlloc(ElementCacheMap::simd_size_double, n_vectors * vec_size * sizeof(double) + n_subsets * sizeof(uint));
	X_ = arena_alloc_->create_array<double>(3*vec_size);
	x_ = X_ + 0;
	y_ = X_ + vec_size;
	z_ = X_ + 2*vec_size;
	if (use_depth_var) d_ = arena_alloc_->create_array<double>(vec_size);
	res_ = arena_alloc_->create_array<double>(vec_size);
	subsets_ = arena_alloc_->create_array<uint>(n_subsets);
    for(unsigned int row=0; row < this->value_.n_rows(); row++)
        for(unsigned int col=0; col < this->value_.n_cols(); col++) {
            // set expression and data to BParser
            unsigned int i_p = row*this->value_.n_cols()+col;
            //b_parser_[i_p].parse(formula_matrix_.at(row,col));
            std::string expr = formula_matrix_.at(row,col); // Need replace some operations to make them compatible with BParser.
                                                            // It will be solved by conversion script after remove fparser, but
                                                            // we mix using of BParser and fparser and need this solution now.
            boost::replace_all(expr, "^", "**"); // power function
            boost::replace_all(expr, "max(", "maximum("); // max function
            boost::replace_all(expr, "min(", "minimum("); // min function
            boost::replace_all(expr, "Pi", "pi"); // Math.pi
            boost::replace_all(expr, "E", "e"); // Math.e
            {  // ternary operator
                std::string pref("if(");
                auto res = std::mismatch(pref.begin(), pref.end(), expr.begin());
                if ( (res.first == pref.end()) && (expr.back() == ')') ) {
                    std::string subexpr = expr.substr(3, expr.size()-4);
                    std::string delimiter = ",";
                    std::string cond = subexpr.substr(0, subexpr.find(delimiter));
                    subexpr.erase(0, cond.size()+1);
                    std::string if_case = subexpr.substr(0, subexpr.find(delimiter));
                    std::string else_case = subexpr.substr(if_case.size()+1);
                    expr = "(" + if_case + " if " + cond + " else " + else_case +")";
                }
            }
            b_parser_[i_p].parse( expr );
            b_parser_[i_p].set_variable("x",  {}, x_);
            b_parser_[i_p].set_variable("y",  {}, y_);
            b_parser_[i_p].set_variable("z",  {}, z_);
            if (use_depth_var) b_parser_[i_p].set_variable("d",  {}, d_);
            b_parser_[i_p].set_constant("t",  {}, {this->time_.end()});
            b_parser_[i_p].set_variable("_result_", {}, res_);
            b_parser_[i_p].compile();
        }
    for (uint i=0; i<n_subsets; ++i)
        subsets_[i] = i;
}


template <int spacedim, class Value>
inline arma::vec FieldFormula<spacedim, Value>::eval_depth_var(const Point &p)
{
	if (surface_depth_ && has_depth_var_) {
		// add value of depth
		arma::vec p_depth(spacedim+1);
		p_depth.subvec(0,spacedim-1) = p;
		try {
			p_depth(spacedim) = surface_depth_->compute_distance(p);
		} catch (SurfaceDepth::ExcTooLargeSnapDistance &e) {
			e << SurfaceDepth::EI_FieldTime(this->time_.end());
			e << in_rec_.ei_address();
			throw;
		}
		return p_depth;
	} else {
		return p;
	}
}


template <int spacedim, class Value>
void FieldFormula<spacedim, Value>::set_dependency(FieldSet &field_set) {
    field_set_ = &field_set;
    for(auto field : field_set.field_list) field_set_names_.insert( field->name() );
}


template <int spacedim, class Value>
FieldFormula<spacedim, Value>::~FieldFormula() {}


// Instantiations of FieldFormula
INSTANCE_ALL(FieldFormula)
