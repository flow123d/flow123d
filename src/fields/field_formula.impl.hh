    /*
 * field_formula.impl.hh
 *
 *  Created on: Jan 2, 2013
 *      Author: jb
 */

#ifndef FIELD_FORMULA_IMPL_HH_
#define FIELD_FORMULA_IMPL_HH_


#include "fields/field_formula.hh"
#include "fparser.hh"
#include "input/input_type.hh"
#include <boost/foreach.hpp>

/// Implementation.

namespace it = Input::Type;

FLOW123D_FORCE_LINK_IN_CHILD(field_formula)


template <int spacedim, class Value>
const Input::Type::Record & FieldFormula<spacedim, Value>::get_input_type()
{
	auto input_type = FieldFormula<spacedim, Value>::get_input_type_instance().make_instance().first;
	return static_cast<const it::Record &>( *input_type );
	/*std::vector<it::TypeBase::ParameterPair> param_vec;
	param_vec.push_back( std::make_pair("element_input_type", boost::make_shared<it::String>()) );
	static auto string_val = StringValue::get_input_type();
	static it::Instance inst = it::Instance(string_val, param_vec).close();
	static auto input_type = inst.make_instance().first.get();

    return it::Record("FieldFormula", FieldAlgorithmBase<spacedim,Value>::template_name()+" Field given by runtime interpreted formula.")
            .derive_from(FieldAlgorithmBase<spacedim, Value>::get_input_type())
            .declare_key("value", *input_type, it::Default::obligatory(),
                                        "String, array of strings, or matrix of strings with formulas for individual "
                                        "entries of scalar, vector, or tensor value respectively.\n"
                                        "For vector values, you can use just one string to enter homogeneous vector.\n"
                                        "For square (($N\\times N$))-matrix values, you can use:\n\n"
                                        " - array of strings of size (($N$)) to enter diagonal matrix\n"
                                        " - array of strings of size (($\\frac12N(N+1)$)) to enter symmetric matrix (upper triangle, row by row)\n"
                                        " - just one string to enter (spatially variable) multiple of the unit matrix.\n"
                                        "Formula can contain variables ```x,y,z,t``` and usual operators and functions." )
			.close(); */
}



template <int spacedim, class Value>
const Input::Type::Instance & FieldFormula<spacedim, Value>::get_input_type_instance() {
	std::vector<it::TypeBase::ParameterPair> param_vec;
	param_vec.push_back( std::make_pair("element_input_type", boost::make_shared<it::String>()) );

	static it::Record rec = it::Record("FieldFormula", FieldAlgorithmBase<spacedim,Value>::template_name()+" Field given by runtime interpreted formula.")
		.derive_from(FieldAlgorithmBase<spacedim, Value>::get_input_type())
		.declare_key("value", StringValue::get_input_type(), it::Default::obligatory(),
									"String, array of strings, or matrix of strings with formulas for individual "
									"entries of scalar, vector, or tensor value respectively.\n"
									"For vector values, you can use just one string to enter homogeneous vector.\n"
									"For square (($N\\times N$))-matrix values, you can use:\n\n"
									" - array of strings of size (($N$)) to enter diagonal matrix\n"
									" - array of strings of size (($\\frac12N(N+1)$)) to enter symmetric matrix (upper triangle, row by row)\n"
									" - just one string to enter (spatially variable) multiple of the unit matrix.\n"
									"Formula can contain variables ```x,y,z,t``` and usual operators and functions." )
		.close();
	return it::Instance(rec, param_vec).close();
}



template <int spacedim, class Value>
const int FieldFormula<spacedim, Value>::registrar =
		Input::register_class< FieldFormula<spacedim, Value>, unsigned int >("FieldFormula") +
		FieldFormula<spacedim, Value>::get_input_type().size();



template <int spacedim, class Value>
FieldFormula<spacedim, Value>::FieldFormula( unsigned int n_comp)
: FieldAlgorithmBase<spacedim, Value>(n_comp),
  formula_matrix_(this->value_.n_rows(), this->value_.n_cols()),
  formula_matrix_helper_(formula_matrix_), parser_matrix_(this->value_.n_rows())
{
    for(unsigned int row=0; row < this->value_.n_rows(); row++) {
        parser_matrix_[row].resize(this->value_.n_cols());
    }
}



template <int spacedim, class Value>
void FieldFormula<spacedim, Value>::init_from_input(const Input::Record &rec) {
    // read formulas form input
    formula_matrix_helper_.init_from_input( rec.val<typename StringValue::AccessType>("value") );
    value_input_address_ = rec.address_string();
}


template <int spacedim, class Value>
bool FieldFormula<spacedim, Value>::set_time(const TimeStep &time) {


    bool any_parser_changed = false;


    std::string vars = string("x,y,z").substr(0, 2*spacedim-1);
    // update parsers
    for(unsigned int row=0; row < this->value_.n_rows(); row++)
        for(unsigned int col=0; col < this->value_.n_cols(); col++) {
            // get all variable names from the formula
            std::vector<std::string> var_list;

            FunctionParser tmp_parser;
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wunused-variable"
            {
                int err=tmp_parser.ParseAndDeduceVariables(formula_matrix_.at(row,col), var_list);
                ASSERT( err != FunctionParser::FP_NO_ERROR, "ParseAndDeduceVariables error: %s\n", tmp_parser.ErrorMsg() );
            }
#pragma GCC diagnostic pop

            bool time_dependent = false;
            BOOST_FOREACH(std::string &var_name, var_list ) {
                if (var_name == std::string("t") ) time_dependent=true;
                else if (var_name == "x" || var_name == "y" || var_name == "z") continue;
                else
                    xprintf(Warn, "Unknown variable '%s' in the  FieldFormula[%d][%d] == '%s'\n at the input address:\n %s \n",
                            var_name.c_str(), row, col, formula_matrix_.at(row,col).c_str(),
                            value_input_address_.c_str() );
            }
            if (time_dependent) {
                parser_matrix_[row][col].AddConstant("t", time.end());
            }

            // TODO:
            // - possibly add user defined constants and units here ...
            // - optimization; possibly parse only if time_dependent  || formula_matrix[][] has changed ...

            if (time_dependent || this->time_ == TimeStep() ) {
                parser_matrix_[row][col].Parse(formula_matrix_.at(row,col), vars);

                if ( parser_matrix_[row][col].GetParseErrorType() != FunctionParser::FP_NO_ERROR ) {
                    xprintf(UsrErr, "ParserError: %s\n in the FieldFormula[%d][%d] == '%s'\n at the input address:\n %s \n",
                        parser_matrix_[row][col].ErrorMsg(),
                        row,col,formula_matrix_.at(row,col).c_str(),
                        value_input_address_.c_str());
                }

                parser_matrix_[row][col].Optimize();
                any_parser_changed = true;
            }


        }

    this->time_=time;
    return any_parser_changed;
}


/**
 * Returns one value in one given point. ResultType can be used to avoid some costly calculation if the result is trivial.
 */
template <int spacedim, class Value>
typename Value::return_type const & FieldFormula<spacedim, Value>::value(const Point &p, const ElementAccessor<spacedim> &elm)
{
    for(unsigned int row=0; row < this->value_.n_rows(); row++)
        for(unsigned int col=0; col < this->value_.n_cols(); col++) {
            this->value_(row,col) = parser_matrix_[row][col].Eval(p.memptr());
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
    ASSERT_EQUAL( point_list.size(), value_list.size() );
    for(unsigned int i=0; i< point_list.size(); i++) {
        Value envelope(value_list[i]);
        ASSERT( envelope.n_rows()==this->value_.n_rows(),
                "value_list[%d] has wrong number of rows: %d; should match number of components: %d\n",
                i, envelope.n_rows(),this->value_.n_rows());

        for(unsigned int row=0; row < this->value_.n_rows(); row++)
            for(unsigned int col=0; col < this->value_.n_cols(); col++) {
                envelope(row,col) = parser_matrix_[row][col].Eval(point_list[i].memptr());
            }
    }
}


template <int spacedim, class Value>
FieldFormula<spacedim, Value>::~FieldFormula() {
}


#endif /* FIELD_FORMULA_IMPL_HH_ */
