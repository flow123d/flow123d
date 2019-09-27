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
 * @file    point_sets.cc
 * @brief
 * @author  David Flanderka
 */

#include <armadillo>
#include "fields/point_sets.hh"
#include "fields/composed_quadrature.hh"
#include "mesh/side_impl.hh"
#include "mesh/sides.h"


/******************************************************************************
 * Implementation of EvalSubset methods.
 */



/******************************************************************************
 * Explicit instantiation of templates
 */

template class EvalSubset<1>;
template class EvalSubset<2>;
template class EvalSubset<3>;
