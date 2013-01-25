/*
 * field_base.cc
 *
 *  Created on: Dec 4, 2012
 *      Author: jb
 */



#include "fields/field_base_impl.hh"
#include "fields/field_python_impl.hh"
#include "fields/field_constant_impl.hh"
#include "fields/field_formula_impl.hh"
#include "fields/field_interpolated_p0_impl.hh"
#include "fields/field_add_potential_impl.hh"
#include "fields/field_elementwise_impl.hh"



// Implementation of FieldCommon

#define INSTANCE_DIM_DEP_VALUES( field, dim_from, dim_to)                                                               \
template class field<dim_from, FieldValue<dim_to>::VectorFixed >;                       \
template class field<dim_from, FieldValue<dim_to>::TensorFixed >;                       \

//first dimension independent values then dimension dependent
#define INSTANCE_TO_ALL(field, dim_from) \
template class field<dim_from, FieldValue<0>::Enum >;                       \
template class field<dim_from, FieldValue<0>::Integer >;                       \
template class field<dim_from, FieldValue<0>::Scalar >;                       \
template class field<dim_from, FieldValue<0>::Vector >;                         \
\
INSTANCE_DIM_DEP_VALUES( field, dim_from, 2) \
INSTANCE_DIM_DEP_VALUES( field, dim_from, 3) \

#define INSTANCE_ALL(field) \
INSTANCE_TO_ALL(field, 0) \
INSTANCE_TO_ALL( field, 1) \
INSTANCE_TO_ALL( field, 2) \
INSTANCE_TO_ALL( field, 3)


INSTANCE_ALL(FieldBase)
INSTANCE_ALL(Field)
INSTANCE_ALL(BCField)
INSTANCE_ALL(FieldConstant)
INSTANCE_ALL(FieldPython)
INSTANCE_ALL(FieldFormula)
INSTANCE_ALL(FieldElementwise)
//INSTANCE_ALL(FieldInterpolatedP0)

template class FieldAddPotential<3, FieldValue<0>::Scalar >;
template class FieldAddPotential<2, FieldValue<0>::Scalar >;

