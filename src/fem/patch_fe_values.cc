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
 * @file    patch_fe_values.cc
 * @brief   Class FEValues calculates finite element data on the actual
 *          cells such as shape function values, gradients, Jacobian of
 *          the mapping from the reference cell etc.
 * @author  Jan Stebel, David Flanderka
 */

#include "fem/patch_fe_values.hh"
#include "fem/op_function_impl.hh"
#include "fem/mapping_p1.hh"
#include "fem/fe_system.hh"
#include "fem/fe_values_map.hh"



// explicit instantiation
template void PatchFEValues<3>::initialize<0>(Quadrature&);
template void PatchFEValues<3>::initialize<1>(Quadrature&);
template void PatchFEValues<3>::initialize<2>(Quadrature&);
template void PatchFEValues<3>::initialize<3>(Quadrature&);

template class PatchFEValues<3>;
