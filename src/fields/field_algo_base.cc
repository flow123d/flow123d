/*
 * field_algo_base.cc
 *
 *  Created on: Dec 4, 2012
 *      Author: jb
 */






/**
 * @file Instantiation of descendants of the FieldBase<...>
 */
#include "fields/field_algo_base.impl.hh"

#include "fields/field_python.impl.hh"
#include "fields/field_constant.impl.hh"
#include "fields/field_formula.impl.hh"
#include "fields/field_interpolated_p0.impl.hh"
#include "fields/field_add_potential.impl.hh"
#include "fields/field_elementwise.impl.hh"
#include "fields/field_fe.impl.hh"


INSTANCE_ALL(FieldAlgorithmBase)
INSTANCE_ALL(FieldConstant)
INSTANCE_ALL(FieldPython)
INSTANCE_ALL(FieldFormula)
INSTANCE_ALL(FieldElementwise)
INSTANCE_ALL(FieldInterpolatedP0)
INSTANCE_ALL(FieldFE)

template class FieldAddPotential<3, FieldValue<0>::Scalar >;
//template class FieldAddPotential<2, FieldValue<0>::Scalar >;


