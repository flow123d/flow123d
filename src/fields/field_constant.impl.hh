/*
 * field_constant.impl.hh
 *
 *  Created on: Dec 15, 2012
 *      Author: jb
 */

#ifndef FIELD_CONSTANT_IMPL_HH_
#define FIELD_CONSTANT_IMPL_HH_

#include "fields/field_constant.hh"
#include "input/input_type.hh"


/// Implementation.

namespace it = Input::Type;


template <int spacedim, class Value>
const Input::Type::Record & FieldConstant<spacedim, Value>::get_input_type(
        Input::Type::AbstractRecord &a_type, const typename Value::ElementInputType *eit
        )
{
    return it::Record("FieldConstant", FieldAlgorithmBase<spacedim,Value>::template_name()+" Field constant in space.")
        .derive_from(a_type)
        .declare_key("value", Value::get_input_type(eit), it::Default::obligatory(),
                                    "Value of the constant field.\n"
                                    "For vector values, you can use scalar value to enter constant vector.\n"
                                    "For square (($N\\\\times N$))-matrix values, you can use:\n"
                                    " - vector of size (($N$)) to enter diagonal matrix\n\n"
                                    " - vector of size (($\\\\frac12N(N+1)$)) to enter symmetric matrix (upper triangle, row by row)\n"
                                    " - scalar to enter multiple of the unit matrix." )
        .allow_auto_conversion("value")
		.close();
}

template <int spacedim, class Value>
const int FieldConstant<spacedim, Value>::registrar =
		Input::register_class< FieldConstant<spacedim, Value>, unsigned int >("FieldConstant");


template <int spacedim, class Value>
FieldConstant<spacedim, Value>::FieldConstant( unsigned int n_comp)
: FieldAlgorithmBase<spacedim, Value>(n_comp)
{}


template <int spacedim, class Value>
FieldConstant<spacedim, Value> &FieldConstant<spacedim, Value>::set_value(const typename Value::return_type &val)
{
    this->r_value_ = val;

    return *this;
}


template <int spacedim, class Value>
void FieldConstant<spacedim, Value>::init_from_input(const Input::Record &rec) {
    this->value_.init_from_input( rec.val<typename Value::AccessType>("value") );
}



/**
 * Returns one value in one given point. ResultType can be used to avoid some costly calculation if the result is trivial.
 */
template <int spacedim, class Value>
typename Value::return_type const & FieldConstant<spacedim, Value>::value(const Point &p, const ElementAccessor<spacedim> &elm)
{
    return this->r_value_;
}



/**
 * Returns std::vector of scalar values in several points at once.
 */
template <int spacedim, class Value>
void FieldConstant<spacedim, Value>::value_list (const std::vector< Point >  &point_list, const ElementAccessor<spacedim> &elm,
                   std::vector<typename Value::return_type>  &value_list)
{
    ASSERT_EQUAL( point_list.size(), value_list.size() );

    for(unsigned int i=0; i< point_list.size(); i++) {
        ASSERT( Value(value_list[i]).n_rows()==this->value_.n_rows(),
                "value_list[%d] has wrong number of rows: %d; should match number of components: %d\n",
                i, Value(value_list[i]).n_rows(),this->value_.n_rows());


        value_list[i]=this->r_value_;
    }
}



template <int spacedim, class Value>
FieldConstant<spacedim, Value>::~FieldConstant() {
}



#endif /* FIELD_CONSTANT_IMPL_HH_ */
