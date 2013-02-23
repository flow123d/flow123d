    /*
 * field_formula_impl.hh
 *
 *  Created on: Jan 2, 2013
 *      Author: jb
 */

#ifndef FIELD_FORMULA_IMPL_HH_
#define FIELD_FORMULA_IMPL_HH_


#include "fields/field_formula.hh"
#include "fparser.hh"
#include "input/input_type.hh"

/// Implementation.

namespace it = Input::Type;

template <int spacedim, class Value>
it::Record FieldFormula<spacedim, Value>::input_type = get_input_type(FieldBase<spacedim,Value>::input_type, NULL);

template <int spacedim, class Value>
Input::Type::Record FieldFormula<spacedim, Value>::get_input_type(
        Input::Type::AbstractRecord &a_type, typename Value::ElementInputType *eit
        )
{
    it::Record type
            = it::Record("FieldFormula", FieldBase<spacedim,Value>::template_name()+" Field given by runtime interpreted formula.")
            .derive_from(a_type)
            .declare_key("value", StringValue::get_input_type(NULL), it::Default::obligatory(),
                                        "String, array of strings, or matrix of strings with formulas for individual "
                                        "entries of scalar, vector, or tensor value respectively.\n"
                                        "For vector values, you can use just one string to enter homogeneous vector.\n"
                                        "For square NxN-matrix values, you can use:\n"
                                        "* array of strings of size N to enter diagonal matrix\n"
                                        "* array of strings of size (N+1)*N/2 to enter symmetric matrix (upper triangle, row by row)\n"
                                        "* just one string to enter (spatially variable) multiple of the unit matrix.\n"
                                        "Formula can contain variables x,y,z,t and usual operators and functions." );

    return type;
}




template <int spacedim, class Value>
FieldFormula<spacedim, Value>::FieldFormula( unsigned int n_comp)
: FieldBase<spacedim, Value>(n_comp),
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

    set_time(this->time_);
}


template <int spacedim, class Value>
bool FieldFormula<spacedim, Value>::set_time(double time) {
    this->time_=time;

    std::string vars = string("x,y,z").substr(0, 2*spacedim-1);
    // update parsers
    for(unsigned int row=0; row < this->value_.n_rows(); row++)
        for(unsigned int col=0; col < this->value_.n_cols(); col++) {

            parser_matrix_[row][col].AddConstant("t", this->time_);
            // possibly add other constants
            parser_matrix_[row][col].Parse(formula_matrix_.at(row,col), vars);

            if ( parser_matrix_[row][col].GetParseErrorType() != FunctionParser::FP_NO_ERROR ) {
                xprintf(UsrErr, "FieldFormula at(%d,%d) error while parsing expression:\n '%s'\nParser error message: %s\n",
                        row,col,formula_matrix_.at(row,col).c_str(),parser_matrix_[row][col].ErrorMsg() );
            }

            parser_matrix_[row][col].Optimize();
        }
    return true; // TODO: check if the fromula contains 't' variable, and return true only in that case
}


/**
 * Returns one value in one given point. ResultType can be used to avoid some costly calculation if the result is trivial.
 */
template <int spacedim, class Value>
typename Value::return_type const & FieldFormula<spacedim, Value>::value(const Point<spacedim> &p, const ElementAccessor<spacedim> &elm)
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
void FieldFormula<spacedim, Value>::value_list (const std::vector< Point<spacedim> >  &point_list, const ElementAccessor<spacedim> &elm,
                   std::vector<typename Value::return_type>  &value_list)
{
    ASSERT_SIZES( point_list.size(), value_list.size() );
    for(unsigned int i=0; i< point_list.size(); i++) {
        Value envelope(value_list[i]);

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
