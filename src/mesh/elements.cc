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
 * @file    elements.cc
 * @ingroup mesh
 * @brief   Various element oriented stuff, should be restricted to purely geometric functions
 */

#include <vector>
#include <string>

#include "system/system.hh"
#include "mesh/accessors.hh"
#include "elements.h"
#include "mesh/mesh.h"
#include "mesh/ref_element.hh"

// following deps. should be removed
//#include "materials.hh"
#include "mesh/accessors.hh"
#include "la/distribution.hh"



Element::Element()
: boundary_idx_(NULL),
  neigh_vb(NULL),
  pid_(0),
  n_neighs_vb_(0),
  dim_(0)
{
}


Element::Element(unsigned int dim, RegionIdx reg)
{
    init(dim, reg);
}



void Element::init(unsigned int dim, RegionIdx reg) {
    pid_=0;
    n_neighs_vb_=0;
    neigh_vb=NULL;
    dim_=dim;
    region_idx_=reg;

    edge_idx_.resize( n_sides() );
    permutation_idx_.resize( n_sides() );
    boundary_idx_ = NULL;

    for (unsigned int si=0; si<this->n_sides(); si++) {
        edge_idx_[ si ]=Mesh::undef_idx;
        permutation_idx_[si] = Mesh::undef_idx;
    }
}


Element::~Element() {
    // Can not make deallocation here since then resize of
    // vectors of elements deallocates what should be keeped.
}


/**
 * Count element sides of the space dimension @p side_dim.
 */

/* If we use this method, it will be moved to mesh accessor class.
unsigned int Element::n_sides_by_dim(unsigned int side_dim)
{
    if (side_dim == dim()) return 1;

    unsigned int n = 0;
    for (unsigned int i=0; i<n_sides(); i++)
        if (side(i)->dim() == side_dim) n++;
    return n;
}*/



//-----------------------------------------------------------------------------
// vim: set cindent:
