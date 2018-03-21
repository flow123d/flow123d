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
#include <boost/foreach.hpp>

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
			.declare_key("surface_region", it::String(), it::Default("\".BOUNDARY\""),
										"The name of region set considered as the surface. You have to set surface region if you "
										"want to use formula variable ```d```.")
			.declare_key("max_depth", it::Double(0.0), it::Default("1e06"),
										"Maximal value of surface depth. If no intersection is found or computed depth is greater, "
										"we use this value.")
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
  formula_matrix_(this->value_.n_rows(), this->value_.n_cols())
{
    parser_matrix_.resize(this->value_.n_rows());
    for(unsigned int row=0; row < this->value_.n_rows(); row++) {
        parser_matrix_[row].resize(this->value_.n_cols());
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


    std::string vars = string("x,y,z").substr(0, 2*spacedim-1) + string(",d");
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
                OLD_ASSERT( err != FunctionParser::FP_NO_ERROR, "ParseAndDeduceVariables error: %s\n", tmp_parser.ErrorMsg() );
            }
#pragma GCC diagnostic pop

            bool time_dependent = false;
            BOOST_FOREACH(std::string &var_name, var_list ) {
                if (var_name == std::string("t") ) time_dependent=true;
                else if (var_name == std::string("d") ) has_depth_var_=true;
                else if (var_name == "x" || var_name == "y" || var_name == "z") continue;
                else
                	WarningOut().fmt("Unknown variable '{}' in the  FieldFormula[{}][{}] == '{}'\n at the input address:\n {} \n",
                            var_name, row, col, formula_matrix_.at(row,col), value_input_address );
            }

            // Seems that we can not just add 't' constant to tmp_parser, since it was already Parsed.
            parser_matrix_[row][col].AddConstant("Pi", 3.14159265358979323846);
            parser_matrix_[row][col].AddConstant("E", 2.71828182845904523536);
            if (time_dependent) {
                parser_matrix_[row][col].AddConstant("t", time.end());
            }

            // TODO:
            // - possibly add user defined constants and units here ...
            // - optimization; possibly parse only if time_dependent  || formula_matrix[][] has changed ...
            //parser_matrix_[row][col] = tmp_parser;
            if (time_dependent || this->time_ == TimeStep() ) {
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

    this->time_=time;
    return any_parser_changed;
}


template <int spacedim, class Value>
void FieldFormula<spacedim, Value>::set_mesh(const Mesh *mesh, bool boundary_domain) {
    // create SurfaceDepth object on surface region
	std::string surface_region = in_rec_.val<std::string>("surface_region");
	surface_depth_ = std::make_shared<SurfaceDepth>(mesh, surface_region, in_rec_.val<std::string>("surface_direction"));
	max_depth_ = in_rec_.val<double>("max_depth");
}


/**
 * Returns one value in one given point. ResultType can be used to avoid some costly calculation if the result is trivial.
 */
template <int spacedim, class Value>
typename Value::return_type const & FieldFormula<spacedim, Value>::value(const Point &p, const ElementAccessor<spacedim> &elm)
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
void FieldFormula<spacedim, Value>::value_list (const std::vector< Point >  &point_list, const ElementAccessor<spacedim> &elm,
                   std::vector<typename Value::return_type>  &value_list)
{
	OLD_ASSERT_EQUAL( point_list.size(), value_list.size() );
    for(unsigned int i=0; i< point_list.size(); i++) {
        Value envelope(value_list[i]);
        OLD_ASSERT( envelope.n_rows()==this->value_.n_rows(),
                "value_list[%d] has wrong number of rows: %d; should match number of components: %d\n",
                i, envelope.n_rows(),this->value_.n_rows());
        auto p_depth = this->eval_depth_var(point_list[i]);

        for(unsigned int row=0; row < this->value_.n_rows(); row++)
            for(unsigned int col=0; col < this->value_.n_cols(); col++) {
                envelope(row,col) = this->unit_conversion_coefficient_ * parser_matrix_[row][col].Eval(p_depth.memptr());
            }
    }
}


template <int spacedim, class Value>
arma::vec FieldFormula<spacedim, Value>::eval_depth_var(const Point &p)
{
	arma::vec p_depth(spacedim+1);
	for (unsigned int i=0; i<spacedim; i++) p_depth(i) = p(i);
	if (has_depth_var_) {
		// add value of depth
		p_depth(spacedim) = std::min( surface_depth_->compute_distance(p), max_depth_ );
	} else {
		p_depth(spacedim) = 0;
	}
	return p_depth;
}


template <int spacedim, class Value>
FieldFormula<spacedim, Value>::~FieldFormula() {
}


// Instantiations of FieldFormula
INSTANCE_ALL(FieldFormula)
