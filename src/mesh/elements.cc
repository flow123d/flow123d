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
#include "mesh/side_impl.hh"
#include "elements.h"
#include "mesh/mesh.h"
#include "mesh/ref_element.hh"

// following deps. should be removed
#include "mesh/boundaries.h"
//#include "materials.hh"
#include "mesh/accessors.hh"
#include "la/distribution.hh"



Element::Element()
: node(NULL),
  boundary_idx_(NULL),
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

    node = new Node * [ n_nodes()];
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
 * SET THE "METRICS" FIELD IN STRUCT ELEMENT
 */
double Element::measure() const {
    switch (dim()) {
        case 0:
            return 1.0;
            break;
        case 1:
            return arma::norm(*(node[ 1 ]) - *(node[ 0 ]) , 2);
            break;
        case 2:
            return
                arma::norm(
                    arma::cross(*(node[1]) - *(node[0]), *(node[2]) - *(node[0])),
                    2
                ) / 2.0 ;
            break;
        case 3:
            return fabs(
                arma::dot(
                    arma::cross(*node[1] - *node[0], *node[2] - *node[0]),
                    *node[3] - *node[0] )
                ) / 6.0;
            break;
    }
    return 1.0;
}

double Element::tetrahedron_jacobian() const
{
    OLD_ASSERT(dim_ == 3, "Cannot provide Jacobian for dimension other than 3.");
    return arma::dot( arma::cross(*node[1] - *node[0], 
                                  *node[2] - *node[0]),
                    *node[3] - *node[0]
                    );
}


/**
 * SET THE "CENTRE[]" FIELD IN STRUCT ELEMENT
 */

arma::vec3 Element::centre() const {
    arma::vec3 centre;
    centre.zeros();

    for (unsigned int li=0; li<this->n_nodes(); li++) {
        centre += node[ li ]->point();
    }
    centre /= (double) n_nodes();
    return centre;
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


double Element::quality_measure_smooth(SideIter side) const {
    if (dim_==3) {
        double sum_faces=0;
        double face[4];
        for(unsigned int i=0; i<4; i++, ++side) sum_faces+=( face[i]=side->measure());

        double sum_pairs=0;
        for(unsigned int i=0;i<3;i++)
            for(unsigned int j=i+1;j<4;j++) {
                unsigned int i_line = RefElement<3>::line_between_faces(i,j);
                arma::vec line = *node[RefElement<3>::interact(Interaction<0,1>(i_line))[1]] - *node[RefElement<3>::interact(Interaction<0,1>(i_line))[0]];
                sum_pairs += face[i]*face[j]*arma::dot(line, line);
            }
        double regular = (2.0*sqrt(2.0/3.0)/9.0); // regular tetrahedron
        return fabs( measure()
                * pow( sum_faces/sum_pairs, 3.0/4.0))/ regular;

    }
    if (dim_==2) {
        return fabs(
                measure()/
                pow(
                         arma::norm(*node[1] - *node[0], 2)
                        *arma::norm(*node[2] - *node[1], 2)
                        *arma::norm(*node[0] - *node[2], 2)
                        , 2.0/3.0)
               ) / ( sqrt(3.0) / 4.0 ); // regular triangle
    }
    return 1.0;
}


void Element::get_bounding_box(BoundingBox &bounding_box) const
{
	bounding_box = BoundingBox( this->node[0]->point() );

	for (unsigned int i=1; i<n_nodes(); i++)
		bounding_box.expand( this->node[i]->point() );
}

/*BoundingBox &Element::get_bounding_box_fast(BoundingBox &bounding_box) const
{
    ASSERT_GT(mesh_->element_box_.size(), 0);
    return mesh_->element_box_[mesh_->element.index(this)];
}*/


/*unsigned int Element::get_proc() const
{
  return mesh_->get_el_ds()->get_proc(mesh_->get_row_4_el()[index()]);
}*/
    

//-----------------------------------------------------------------------------
// vim: set cindent:
