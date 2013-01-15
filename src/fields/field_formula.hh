/*
 * field_formula.hh
 *
 *  Created on: Jan 2, 2013
 *      Author: jb
 */

#ifndef FIELD_FORMULA_HH_
#define FIELD_FORMULA_HH_


#include "system/system.hh"
#include "fields/field_base.hh"
#include "mesh/point.hh"

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
class FieldFormula : public FieldBase<spacedim, Value>
{
public:

    FieldFormula(const double init_time=0.0, unsigned int n_comp=0);


    static Input::Type::Record input_type;

    static Input::Type::Record get_input_type(Input::Type::AbstractRecord &a_type, typename Value::ElementInputType *eit);

    virtual void init_from_input( Input::Record rec);

    virtual void set_time(double time);

    /**
     * Returns one value in one given point. ResultType can be used to avoid some costly calculation if the result is trivial.
     */
    virtual typename Value::return_type &value(const Point<spacedim> &p, ElementAccessor<spacedim> &elm);

    /**
     * Returns std::vector of scalar values in several points at once.
     */
    virtual void value_list (const std::vector< Point<spacedim> >  &point_list, ElementAccessor<spacedim> &elm,
                       std::vector<typename Value::return_type>  &value_list);


    virtual ~FieldFormula();

private:
    typedef FieldValue_<Value::NRows_, Value::NCols_, std::string> StringValue;

    typename StringValue::return_type formula_matrix_;
    StringValue formula_matrix_helper_;
    std::vector< std::vector<FunctionParser> > parser_matrix_;

};



#endif /* FIELD_FORMULA_HH_ */
