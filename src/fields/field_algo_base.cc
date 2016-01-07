/*!
 *
﻿ * Copyright (C) 2015 Technical University of Liberec.  All rights reserved.
 * 
 * This program is free software; you can redistribute it and/or modify it under
 * the terms of the GNU General Public License version 3 as published by the
 * Free Software Foundation. (http://www.gnu.org/licenses/gpl-3.0.en.html)
 * 
 * This program is distributed in the hope that it will be useful, but WITHOUT
 * ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 * FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.
 *
 * 
 * @file    field_algo_base.cc
 * @brief   
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


