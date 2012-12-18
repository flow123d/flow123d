/*
 * field_constant_impl.hh
 *
 *  Created on: Dec 15, 2012
 *      Author: jb
 */

#include "fields/field_all.hh"

#ifndef FIELD_CONSTANT_IMPL_HH_
#define FIELD_CONSTANT_IMPL_HH_

#include "fields/field_constant.hh"


//#include <boost/type_traits.hpp>

/// Implementation.

template <int spacedim, class Value>
FieldConstant<spacedim, Value>::FieldConstant(const double init_time, unsigned int n_comp)
: FieldBase<spacedim, Value>(init_time, n_comp)
{}



template <int spacedim, class Value>
Input::Type::Record &FieldConstant<spacedim, Value>::get_input_type() {
    using namespace  Input::Type;

    static Record rec("FieldConstant", FieldBase<spacedim,Value>::template_name()+" Field constant in space.");

    if (! rec.is_finished()) {
        rec.derive_from(FieldBase<spacedim, Value>::get_input_type());
        rec.declare_key("value", Value::get_input_type(), Default("0.0"),
                                "Value of the constant field.\n"
                                "For vector values, you can use scalar value to enter constatn vector.\n"
                                "For square NxN-matrix values, you can use:\n"
                                "* vector of size N to enter diagonal matrix\n"
                                "* vector of size (N+1)*N/2 to enter symmetric matrix (upper triangle, row by row)\n"
                                "* scalar to enter multiple of the unit matrix." );
        rec.finish();
    }
    return rec;
}



template <int spacedim, class Value>
void FieldConstant<spacedim, Value>::init_from_input( Input::Record rec) {
    this->value_.init_from_input( rec.val<typename Value::InputType>("value") );
}







/**
* Returns one vector value in one given point.
*/
template <int spacedim, class Value>
FieldResult FieldConstant<spacedim, Value>::value(const Point<spacedim> &p, ElementAccessor<spacedim> &elm, typename Value::return_type &value)
{
    value=this->value_;
    return result_other;
}



/*
template <int spacedim, class Value>
void FieldConstant<spacedim, Value>::value_list (const std::vector< Point<spacedim> >  &point_list, ElementAccessor<spacedim> &elm,
                   std::vector<Value>  &value_list, std::vector<FieldResult> &result_list)
{
#ifdef HAVE_PYTHON
    ASSERT_SIZES( point_list.size(), value_list.size() );
    for(unsigned int i=0; i< point_list.size(); i++)
        result_list[i] = value(point_list[i], component);
#endif // HAVE_PYTHON
}
*/



template <int spacedim, class Value>
FieldConstant<spacedim, Value>::~FieldConstant() {
}



#endif /* FIELD_CONSTANT_IMPL_HH_ */
