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
 * @file    bc_multi_field.hh
 * @brief   
 */

#ifndef BC_MULTI_FIELD_HH_
#define BC_MULTI_FIELD_HH_


#include "multi_field.hh"


/**
 * Same as MultiField<...> but for boundary regions.
 */
template<int spacedim, class Value>
class BCMultiField : public MultiField<spacedim, Value> {
public:
    BCMultiField() : MultiField<spacedim,Value>(true) {}
};

#endif /* BC_MULTI_FIELD_HH_ */
