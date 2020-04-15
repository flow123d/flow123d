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
 * @file    discrete_space.cc
 * @brief   Implementation of class which provides the finite element for every mesh cell.
 * @author  Jan Stebel
 */

#include "fem/discrete_space.hh"
#include "fem/finite_element.hh"
#include "mesh/accessors.hh"

MixedPtr<FiniteElement> EqualOrderDiscreteSpace::fe(FMT_UNUSED const ElementAccessor<3> &cell) const { return fe_; }





