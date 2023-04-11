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
#include "fields/field_instances.hh"	// for instantiation macros

#include "fields/field_python.impl.hh"


FLOW123D_FORCE_LINK_IN_PARENT(field_constant)
FLOW123D_FORCE_LINK_IN_PARENT(field_formula)
FLOW123D_FORCE_LINK_IN_PARENT(field_python)
FLOW123D_FORCE_LINK_IN_PARENT(field_time_function)
FLOW123D_FORCE_LINK_IN_PARENT(field_fe)


INSTANCE_ALL(FieldAlgorithmBase)
INSTANCE_ALL(FieldPython)

//template class FieldAddPotential<3, FieldValue<0>::Scalar >;
//template class FieldAddPotential<2, FieldValue<0>::Scalar >;

// temporary solution for computing more fields at once in python
template class FieldPython<3, FieldValue<0>::Vector >;
