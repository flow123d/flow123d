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
 * @file    field_values.cc
 * @brief   
 */

#include "fields/field_values.hh"

namespace internal {
// Helper functions to get scalar type name
std::string type_name_(double) { return "Real"; }
std::string type_name_(int)  { return "Int"; }
std::string type_name_(FieldEnum) { return "Enum"; }


} // namespace internal


template class FieldValue_<1,1,FieldEnum>;
template class FieldValue_<1,1,int>;
template class FieldValue_<0,1,int>;
template class FieldValue_<1,1,double>;
template class FieldValue_<0,1,double>;

template class FieldValue_<2,1,double>;
template class FieldValue_<3,1,double>;

template class FieldValue_<2,2,double>;
template class FieldValue_<3,3,double>;

template class FieldValue<1>;
template class FieldValue<2>;
template class FieldValue<3>;
