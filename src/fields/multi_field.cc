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
 * @file    multi_field.cc
 * @brief   
 */

#include "fields/field_algo_base.impl.hh"	// for instantiation macros

#include "fields/multi_field.hh"
#include "fields/multi_field.impl.hh"



/****************************************************************************
 *  Instances of field templates
 */



template class MultiField<3, FieldValue<0>::Scalar >;
template class MultiField<3, FieldValue<3>::VectorFixed >;
template class MultiField<3, FieldValue<3>::TensorFixed >;
template class MultiField<3, FieldValue<0>::Enum >;
template class MultiField<3, FieldValue<0>::Integer >;
//template class MultiField<2, FieldValue<0>::Scalar >;
//template class MultiField<2, FieldValue<0>::Enum >;
//template class MultiField<2, FieldValue<0>::Integer >;
