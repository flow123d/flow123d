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
 * @file    bc_field.hh
 * @brief   
 */

#ifndef BC_FIELD_HH_
#define BC_FIELD_HH_


#include "field.hh"


/**
 * Same as Field<...> but for boundary regions.
 *
 * Definition of BCField must be in separate file.
 * In other case source file field.cc is too big and compiler can throw compile error.
 */
template<int spacedim, class Value>
class BCField : public Field<spacedim, Value> {
public:
    BCField() : Field<spacedim,Value>("anonymous_bc", true) {}
};

#endif /* BC_FIELD_HH_ */
