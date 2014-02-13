/*
 * field_base.cc
 *
 *  Created on: Dec 4, 2012
 *      Author: jb
 */




/**
 * @file Instantiation of descendants of the FieldBase<...>
 */
#include "fields/field_base_impl.hh"

#include "fields/field_python_impl.hh"
#include "fields/field_constant_impl.hh"
#include "fields/field_formula_impl.hh"
#include "fields/field_interpolated_p0_impl.hh"
#include "fields/field_add_potential_impl.hh"
#include "fields/field_elementwise_impl.hh"


INSTANCE_ALL(FieldBase)
INSTANCE_ALL(FieldConstant)
INSTANCE_ALL(FieldPython)
INSTANCE_ALL(FieldFormula)
INSTANCE_ALL(FieldElementwise)
INSTANCE_ALL(FieldInterpolatedP0)

template class FieldAddPotential<3, FieldValue<0>::Scalar >;
//template class FieldAddPotential<2, FieldValue<0>::Scalar >;


