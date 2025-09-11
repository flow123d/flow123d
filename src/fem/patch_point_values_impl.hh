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
 * @file    patch_point_values_impl.hh
 * @brief   Store finite element data on the actual patch
 *          such as shape function values, gradients, Jacobian
 *          of the mapping from the reference cell etc.
 * @author  David Flanderka
 */

#ifndef PATCH_POINT_VALUES_HH_IMPL_
#define PATCH_POINT_VALUES_HH_IMPL_

#include "fem/patch_point_values.hh"


namespace Op {
    class BulkDomain;
    class SideDomain;
}

// Template specialized methods

template<>
template<>
NodeAccessor<3> PatchPointValues<3>::node<Op::BulkDomain>(unsigned int i_elm, unsigned int i_n) {
    return (*elem_dim_list_)[i_elm].node(i_n);
}

template<>
template<>
NodeAccessor<3> PatchPointValues<3>::node<Op::SideDomain>(unsigned int i_elm, unsigned int i_n) {
    return side_list_[i_elm].node(i_n);
}

template<>
template<>
unsigned int PatchPointValues<3>::n_mesh_entities<Op::BulkDomain>() {
    return elem_dim_list_->size();
}

template<>
template<>
unsigned int PatchPointValues<3>::n_mesh_entities<Op::SideDomain>() {
    return side_list_.size();
}



#endif /* PATCH_POINT_VALUES_HH_IMPL_ */
