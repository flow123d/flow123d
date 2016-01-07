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
#include "mesh/mesh.h"
#include "mesh/ref_element.hh"
#include "element_impls.hh"

// following deps. should be removed
#include "mesh/boundaries.h"
//#include "materials.hh"
#include "mesh/accessors.hh"



Element::Element()
:  pid(0),

  node(NULL),

//  material(NULL),
  edge_idx_(NULL),
  boundary_idx_(NULL),
  permutation_idx_(NULL),

  n_neighs_vb(0),
  neigh_vb(NULL),

  dim_(0)

{
}


Element::Element(unsigned int dim, Mesh *mesh_in, RegionIdx reg)
{
    init(dim, mesh_in, reg);
}



void Element::init(unsigned int dim, Mesh *mesh_in, RegionIdx reg) {
    pid=0;
    n_neighs_vb=0;
    neigh_vb=NULL;
    dim_=dim;
    mesh_=mesh_in;
    region_idx_=reg;

    node = new Node * [ n_nodes()];
    edge_idx_ = new unsigned int [ n_sides()];
    boundary_idx_ = NULL;
    permutation_idx_ = new unsigned int[n_sides()];

    FOR_ELEMENT_SIDES(this, si) {
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
    ASSERT(dim_ == 3, "Cannot provide Jacobian for dimension other than 3.");
    return arma::dot( arma::cross(*node[1] - *node[0], 
                                  *node[2] - *node[0]),
                    *node[3] - *node[0]
                    );
}


/**
 * SET THE "CENTRE[]" FIELD IN STRUCT ELEMENT
 */

arma::vec3 Element::centre() const {
    unsigned int li;

    arma::vec3 centre;
    centre.zeros();

    FOR_ELEMENT_NODES(this, li) {
        centre += node[ li ]->point();
    }
    centre /= (double) n_nodes();
    return centre;
}

/**
 * Count element sides of the space dimension @p side_dim.
 */

unsigned int Element::n_sides_by_dim(unsigned int side_dim)
{
    if (side_dim == dim()) return 1;

    unsigned int n = 0;
    for (unsigned int i=0; i<n_sides(); i++)
        if (side(i)->dim() == side_dim) n++;
    return n;
}


ElementAccessor< 3 > Element::element_accessor() const
{
  return mesh_->element_accessor( mesh_->element.index(this) );
}



Region Element::region() const {
    return Region( region_idx_, mesh_->region_db());
}


unsigned int Element::id() const {
	return mesh_->element.get_id(this);
}

double Element::quality_measure_smooth() {
    if (dim_==3) {
        double sum_faces=0;
        double face[4];
        for(unsigned int i=0;i<4;i++) sum_faces+=( face[i]=side(i)->measure());

        double sum_pairs=0;
        for(unsigned int i=0;i<3;i++)
            for(unsigned int j=i+1;j<4;j++) {
                unsigned int i_line = RefElement<3>::line_between_faces(i,j);
//                 arma::vec line = *node[RefElement<3>::line_nodes[i_line][1]] - *node[RefElement<3>::line_nodes[i_line][0]];
                arma::vec line = *node[RefElement<3>::interact<0,1>(i_line)[1]] - *node[RefElement<3>::interact<0,1>(i_line)[0]];
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

//-----------------------------------------------------------------------------
// vim: set cindent:
