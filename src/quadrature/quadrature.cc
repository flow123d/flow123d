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
 * @file    quadrature.hh
 * @brief   Basic definitions of numerical quadrature rules.
 * @author  Jan Stebel
 */

#include "quadrature/quadrature.hh"



Quadrature& Quadrature::operator=(const Quadrature &q)
{
    ASSERT_DBG( dim_ == q.dim_ );
    quadrature_points = q.quadrature_points;
    weights = q.weights;
    return *this;
}
    

Quadrature::Quadrature(unsigned int dimension, unsigned int n_q)
: dim_(dimension),
  quadrature_points(dimension, 1,  n_q),
  weights(n_q, 0)
{}


Quadrature::Quadrature(const Quadrature &q) :
        dim_(q.dim_),
        quadrature_points(q.quadrature_points),
        weights(q.weights)
{}


template<unsigned int bulk_dim>
Quadrature Quadrature::make_from_side(unsigned int sid, unsigned int pid) const
{
    ASSERT_DBG( bulk_dim == dim_ + 1 );
    
    // Below we permute point coordinates according to permutation
    // of nodes on side. We just check that these numbers equal.
    ASSERT_DBG( RefElement<bulk_dim>::n_nodes_per_side == bulk_dim );
    
    Quadrature q(dim_+1, size());
    Armor::ArmaVec<double, bulk_dim+1> el_bar_coords, final_bar;
    
    for (unsigned int k=0; k<size(); k++)
    {
        //compute barycentric coordinates on element
        Armor::ArmaVec<double, bulk_dim> p = RefElement<bulk_dim-1>::local_to_bary(point<bulk_dim-1>(k));
        Armor::ArmaVec<double, bulk_dim> pp;
        
        //permute
        for (unsigned int i=0; i<RefElement<bulk_dim>::n_nodes_per_side; i++) {
            pp(RefElement<bulk_dim>::side_permutations[pid][i]) = p(i);
        }
        
        el_bar_coords = RefElement<bulk_dim>::template interpolate<bulk_dim-1>(pp, sid);
        
        //get local coordinates and set
        q.quadrature_points.set(k) = RefElement<bulk_dim>::bary_to_local(el_bar_coords);
        q.weights[k] = weight(k);
    }
    
    return q;
}

// Specialized subquadrature consructor for dim=1.
template<> Quadrature Quadrature::make_from_side<1>(unsigned int sid, unsigned int) const
{
    ASSERT_EQ_DBG(size(), 1);
    Quadrature q(1, 1);
    q.quadrature_points.set(0) = Armor::ArmaVec<double, 1>({ (double)sid });
    q.weight(0) = 1;

    return q;
}

template Quadrature Quadrature::make_from_side<2>(unsigned int sid, unsigned int pid) const;
template Quadrature Quadrature::make_from_side<3>(unsigned int sid, unsigned int pid) const;


