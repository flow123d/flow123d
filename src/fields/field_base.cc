/*
 * field_base.cc
 *
 *  Created on: Dec 4, 2012
 *      Author: jb
 */



#include "fields/field_base_impl.hh"
#include "fields/field_python_impl.hh"
#include "fields/field_constant_impl.hh"
#include "fields/field_interpolated_p0_impl.hh"

/*
#define INSTANCE_FIELD_ALL_TYPES(field)    \
template class field<3, FieldValue<0>::Discrete >;  \
template class field<3, FieldValue<0>::Scalar >;    \
template class field<3, FieldValue<0>::Vector >;        \
template class field<3, FieldValue<3>::VectorFixed >;    \
template class field<3, FieldValue<3>::TensorFixed >;     \
*/


// Implementation of FieldCommon

#define INSTANCE_DIM_DEP_VALUES( field, dim_from, dim_to)                                                               \
template class field<dim_from, FieldValue<dim_to>::VectorFixed >;                       \
template class field<dim_from, FieldValue<dim_to>::TensorFixed >;                       \

//first dimension independent values then dimension dependent
#define INSTANCE_TO_ALL(field, dim_from) \
template class field<dim_from, FieldValue<0>::Discrete >;                       \
template class field<dim_from, FieldValue<0>::Scalar >;                       \
template class field<dim_from, FieldValue<0>::Vector >;                         \
\
INSTANCE_DIM_DEP_VALUES( field, dim_from, 2) \
INSTANCE_DIM_DEP_VALUES( field, dim_from, 3) \

#define INSTANCE_ALL(field) \
INSTANCE_TO_ALL(field, 0) \
INSTANCE_TO_ALL( field, 1) \
INSTANCE_TO_ALL( field, 2) \
INSTANCE_TO_ALL( field, 3) \


INSTANCE_ALL(FieldBase)
INSTANCE_ALL(FieldConstant)
INSTANCE_ALL(FieldPython)
//INSTANCE_ALL(FieldInterpolatedP0)
