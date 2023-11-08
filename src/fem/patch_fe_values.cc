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


template<unsigned int spacedim>
void PatchFEValues<spacedim>::resize(unsigned int new_size) {
    ASSERT_LE(new_size, max_size());
    used_size_ = new_size;

    // TODO maybe reset Element data objects??
}

template<unsigned int spacedim>
void PatchFEValues<spacedim>::allocate_in()
{
//    if (this->update_flags & update_values)
//        shape_values_.resize(this->n_points_, vector<double>(this->n_dofs_*this->n_components_));
//
//    if (this->update_flags & update_gradients)
//        shape_gradients_.resize(this->n_points_, vector<arma::vec::fixed<spacedim> >(this->n_dofs_*this->n_components_));

    this->fv_ = this;
}


template<unsigned int spacedim>
void PatchFEValues<spacedim>::initialize_in(
         FMT_UNUSED Quadrature &q,
		 FMT_UNUSED unsigned int dim)
{
}

/// Explicit initialization
template class PatchFEValues<3>;
