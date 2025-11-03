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
#include "fem/patch_point_values.hh"
#include "fem/op_accessors_impl.hh"
#include "fem/mapping_p1.hh"
#include "fem/fe_system.hh"
#include "fem/fe_values_map.hh"



template<unsigned int spacedim>
template<unsigned int dim>
unsigned int PatchFEValues<spacedim>::n_dofs_high() const {
    ASSERT((dim>=0) && (dim<=2))(dim).error("Dimension must be 0, 1, 2.");
    return fe_[Dim<dim+1>{}]->n_dofs();
}

template<>
template<>
unsigned int PatchFEValues<3>::n_dofs_high<3>() const {
    return fe_[Dim<3>{}]->n_dofs();
}




// explicit instantiation
template unsigned int PatchFEValues<3>::n_dofs_high<1>() const;
template unsigned int PatchFEValues<3>::n_dofs_high<2>() const;
template unsigned int PatchFEValues<3>::n_dofs_high<3>() const;

template class PatchFEValues<3>;
