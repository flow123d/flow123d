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
    

Quadrature::Quadrature(unsigned int dimension, const unsigned int n_q)
: dim_(dimension),
  quadrature_points(n_q, dimension),
  weights(n_q, 0)
{}


Quadrature::Quadrature(const Quadrature &q) :
        dim_(q.dim_),
        quadrature_points(q.quadrature_points),
        weights(q.weights)
{}







template<unsigned int quad_dim>
Quadrature quadrature_from_side(const Quadrature &subq, unsigned int sid, unsigned int pid)
{
    ASSERT_DBG( subq.dim() + 1 == quad_dim );
    
    // Below we permute point coordinates according to permutation
    // of nodes on side. We just check that these numbers equal.
    ASSERT_DBG( RefElement<quad_dim>::n_nodes_per_side == quad_dim );
    
    Quadrature q(subq.dim()+1, subq.size());

//     double lambda;
// 
//     // vectors of barycentric coordinates of quadrature points
//     arma::vec::fixed<dim+1> el_bar_coords;
//     arma::vec::fixed<dim> side_bar_coords;
// 
//     for (unsigned int k=0; k<subq.size(); k++)
//     {
//         const arma::vec::fixed<dim-1> &sub_point = subq.point(k);
//         // Calculate barycentric coordinates on the side of the k-th
//         // quadrature point.
//         el_bar_coords.zeros();
//         lambda = 0;
//         // Apply somewhere permutation of indices!
//         for (unsigned int j=0; j<dim-1; j++)
//         {
//             side_bar_coords(j) = sub_point(j);
//             lambda += sub_point(j);
//         }
//         side_bar_coords(dim-1) = 1.0 - lambda;
// 
//         // transform to element coordinates
//         auto side_nodes = RefElement<dim>::interact(Interaction<0, (dim - 1)>(sid));
//         for (unsigned int i=0; i<dim; i++) {
//             // TODO: use RefElement<>::interpolate to map coordinates from the subelement
//             unsigned int i_node = (side_nodes[RefElement<dim>::side_permutations[pid][i]]+dim)%(dim+1);
//             el_bar_coords(i_node) = side_bar_coords((i+dim-1)%dim);
//         }
//         quadrature_points[k] = el_bar_coords.subvec(0,dim-1);
//         weights[k] = subq.weight(k);
//     }
    
    arma::vec::fixed<quad_dim+1> el_bar_coords, final_bar;
    
    for (unsigned int k=0; k<subq.size(); k++)
    {
        //compute barycentric coordinates on element
        arma::vec::fixed<quad_dim> p = RefElement<quad_dim-1>::local_to_bary(subq.point<quad_dim-1>(k).arma());
        arma::vec::fixed<quad_dim> pp;
        
        //permute
        for (unsigned int i=0; i<RefElement<quad_dim>::n_nodes_per_side; i++) {
            pp(RefElement<quad_dim>::side_permutations[pid][i]) = p(i);
        }
        
        el_bar_coords = RefElement<quad_dim>::template interpolate<quad_dim-1>(pp,sid);
        
        //get local coordinates and set
        q.point<quad_dim>(k) = RefElement<quad_dim>::bary_to_local(el_bar_coords);
        q.weight(k) = subq.weight(k);
    }
    
    return q;
}

// Specialized subquadrature consructor for dim=1.
template<> Quadrature quadrature_from_side<1>(const Quadrature &subq, unsigned int sid, unsigned int pid)
{
    Quadrature q(1, 1);
    q.point<1>(0) = { (double)sid };
    q.weight(0) = 1;
    
    return q;
}

template Quadrature quadrature_from_side<2>(const Quadrature &subq, unsigned int sid, unsigned int pid);
template Quadrature quadrature_from_side<3>(const Quadrature &subq, unsigned int sid, unsigned int pid);
