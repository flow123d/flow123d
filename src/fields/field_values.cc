/*
 * field_values.cc
 *
 *  Created on: Dec 11, 2012
 *      Author: jb
 */


#include "fields/field_values.hh"

template class FieldValue<0>::Discrete;
template class FieldValue<0>::Scalar;
template class FieldValue<0>::Vector;

template class FieldValue<2>::VectorFixed;
template class FieldValue<2>::TensorFixed;

template class FieldValue<3>::VectorFixed;
template class FieldValue<3>::TensorFixed;


template class FieldValue<1>;
template class FieldValue<2>;
template class FieldValue<3>;
