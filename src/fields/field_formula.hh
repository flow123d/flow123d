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

/**
 * Class representing fields given by runtime parsed formulas.
 *
 */
template <int spacedim, class Value>
class FieldFormula : public FieldBase<spacedim, Value>
{
public:

    FieldFormula(const double init_time=0.0, unsigned int n_comp=0);

    static Input::Type::Record input_type;

    virtual void init_from_input( Input::Record rec);

    using FieldBase<spacedim,Value>::value;

    /**
     * Returns one value in one given point. ResultType can be used to avoid some costly calculation if the result is trivial.
     */
    virtual FieldResult value(const Point<spacedim> &p, ElementAccessor<spacedim> &elm, typename Value::return_type &value);

    /**
     * Returns std::vector of scalar values in several points at once.
     */
//    virtual void value_list (const std::vector< Point<spacedim> >  &point_list, ElementAccessor<spacedim> &elm,
//                       std::vector<Value>  &value_list,
//                       std::vector<FieldResult> &result_list);


    virtual ~FieldFormula();

private:

};



#endif /* FIELD_FORMULA_HH_ */
