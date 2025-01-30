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
#include "fem/op_factory.hh"
#include "fem/op_function_impl.hh"
#include "fem/op_accessors_impl.hh"
#include "fem/mapping_p1.hh"
#include "fem/fe_system.hh"
#include "fem/fe_values_map.hh"



template<unsigned int spacedim>
template<unsigned int dim>
BulkValues<dim> PatchFEValues<spacedim>::bulk_values() {
  	ASSERT((dim>0) && (dim<=3))(dim).error("Dimension must be 1, 2 or 3.");
    return BulkValues<dim>(&patch_point_vals_bulk_[dim-1], *this, fe_);
}

template<unsigned int spacedim>
template<unsigned int dim>
SideValues<dim> PatchFEValues<spacedim>::side_values() {
   	ASSERT((dim>0) && (dim<=3))(dim).error("Dimension must be 1, 2 or 3.");
    return SideValues<dim>(&patch_point_vals_side_[dim-1], *this, fe_);
}

template<unsigned int spacedim>
template<unsigned int dim>
JoinValues<dim> PatchFEValues<spacedim>::join_values() {
   	//ASSERT((dim>1) && (dim<=3))(dim).error("Dimension must be 2 or 3.");
    return JoinValues<dim>(&patch_point_vals_bulk_[dim-2], &patch_point_vals_side_[dim-1], *this, fe_);
}




// explicit instantiation
template void PatchFEValues<3>::initialize<0>(Quadrature&);
template void PatchFEValues<3>::initialize<1>(Quadrature&);
template void PatchFEValues<3>::initialize<2>(Quadrature&);
template void PatchFEValues<3>::initialize<3>(Quadrature&);
template BulkValues<1> PatchFEValues<3>::bulk_values<1>();
template BulkValues<2> PatchFEValues<3>::bulk_values<2>();
template BulkValues<3> PatchFEValues<3>::bulk_values<3>();
template SideValues<1> PatchFEValues<3>::side_values<1>();
template SideValues<2> PatchFEValues<3>::side_values<2>();
template SideValues<3> PatchFEValues<3>::side_values<3>();
template JoinValues<1> PatchFEValues<3>::join_values<1>();
template JoinValues<2> PatchFEValues<3>::join_values<2>();
template JoinValues<3> PatchFEValues<3>::join_values<3>();

template class PatchFEValues<3>;
