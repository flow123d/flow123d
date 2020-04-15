/*!
 *
﻿ * Copyright (C) 2019 Technical University of Liberec.  All rights reserved.
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
 * @file    fe_values.cc
 * @brief   Class ElementValues calculates data related to transformation
 *          of reference cell to actual cell (Jacobian, inverse Jacobian,
 *          determinant, point coordinates etc.).
 * @author  Jan Stebel
 */

#include "fem/mapping_p1.hh"
#include "quadrature/quadrature.hh"
#include "fem/element_values.hh"



using namespace arma;
using namespace std;








RefElementData::RefElementData(unsigned int np)
    : n_points(np)
{
    bar_coords.resize(np);
    weights.resize(np);
}


template<unsigned int spacedim>
ElementData<spacedim>::ElementData(unsigned int size,
                         UpdateFlags flags,
                         unsigned int dim)
: dim_(dim),
  JxW_values(flags & update_JxW_values ? size : 0),
  side_JxW_values(flags & update_side_JxW_values ? size : 0),
  jacobians(spacedim, dim, flags & update_jacobians ? size : 0),
  determinants(flags & update_volume_elements ? size : 0),
  inverse_jacobians(dim, spacedim, flags & update_inverse_jacobians ? size : 0),
  points(spacedim, 1, flags & update_quadrature_points ? size : 0),
  normal_vectors(spacedim, 1, flags & update_normal_vectors ? size : 0),
  update_flags(flags)
{}


template<unsigned int spacedim>
void ElementData<spacedim>::print()
{
    if (cell.is_valid() || side.is_valid())
    {
        if (cell.is_valid())
            printf("cell %d dim %d ", cell.elm_idx(), cell.dim());
        else if (side.is_valid())
            printf("cell %d dim %d side %d ", side.elem_idx(), side.dim(), side.side_idx());
        
        printf(" det[");
        for (auto d : determinants) printf("%g ", d); printf("]");
        
        printf(" JxW[");
        for (auto j : JxW_values) printf("%g ", j); printf("]");

        printf(" side_JxW[");
        for (auto j : side_JxW_values) printf("%g ", j); printf("]");
        
        printf(" nv[");
        for (unsigned int i=0; i<normal_vectors.size(); i++)
        {
            auto n = normal_vectors.vec<spacedim>(i);
            printf(" [");
            for (unsigned int c=0; c<spacedim; c++) printf("%g ", n[c]);
            printf("]");
        }
        printf("]\n");
    }
}





template<unsigned int spacedim>
RefElementData *ElementValues<spacedim>::init_ref_data(const Quadrature &q)
{
    ASSERT_DBG( q.dim() == dim_ );
    ASSERT_DBG( q.size() == n_points_ );
    RefElementData *ref_data = new RefElementData(q.size());

    for (unsigned int i=0; i<q.size(); i++)
    {
        switch (q.dim())
        {
            case 1:
                ref_data->bar_coords[i] = RefElement<1>::local_to_bary(q.point<1>(i));
                break;
            case 2:
                ref_data->bar_coords[i] = RefElement<2>::local_to_bary(q.point<2>(i));
                break;
            case 3:
                ref_data->bar_coords[i] = RefElement<3>::local_to_bary(q.point<3>(i));
                break;
        default:
            ASSERT(false)(q.dim()).error("Unsupported dimension.\n");
            break;
        }
        ref_data->weights[i] = q.weight(i);
    }
    
    return ref_data;
}



template<unsigned int spacedim>
UpdateFlags ElementValues<spacedim>::update_each(UpdateFlags flags)
{
    switch (dim_)
    {
        case 0:
            flags = MappingP1<0,spacedim>::update_each(flags);
            break;
        case 1:
            flags = MappingP1<1,spacedim>::update_each(flags);
            break;
        case 2:
            flags = MappingP1<2,spacedim>::update_each(flags);
            break;
        case 3:
            flags = MappingP1<3,spacedim>::update_each(flags);
            break;
        default:
            ASSERT(false)(dim_).error("Unsupported dimension.\n");
            break;
    }
    return flags;
}



template<unsigned int spacedim>
ElementValues<spacedim>::ElementValues(
         Quadrature &_quadrature,
         UpdateFlags _flags,
         unsigned int dim)
: dim_(dim),
  n_points_(_quadrature.size()),
  n_sides_(_quadrature.dim() == dim ? 0 : dim+1),
  n_side_permutations_(_quadrature.dim() +1 == dim ? (dim+1)*(2*dim*dim-5*dim+6)/6 : 0),
  ref_data(nullptr),
  side_ref_data(n_sides_, std::vector<RefElementData*>(n_side_permutations_)),
  data(n_points_, update_each(_flags), dim)
{
    if (dim == 0) return; // avoid unnecessary allocation of dummy 0 dimensional objects
    if ( _quadrature.dim() == dim )
    {
        // precompute element data
        ref_data = init_ref_data(_quadrature);
    }
    else if ( _quadrature.dim() + 1 == dim )
    {
        // precompute side data
        for (unsigned int sid = 0; sid < n_sides_; sid++)
        {
            for (unsigned int pid = 0; pid < n_side_permutations_; pid++)
            {
                Quadrature side_quad(dim);
                // transform the side quadrature points to the cell quadrature points
                switch (dim)
                {
                    case 1:
                        side_quad = _quadrature.make_from_side<1>(sid, pid);
                        break;
                    case 2:
                        side_quad = _quadrature.make_from_side<2>(sid, pid);
                        break;
                    case 3:
                        side_quad = _quadrature.make_from_side<3>(sid, pid);
                        break;
                    default:
                        ASSERT(false)(dim).error("Unsupported dimension.\n");
                        break;
                }
                side_ref_data[sid][pid] = init_ref_data(side_quad);
            }
        }
    }
}


template<unsigned int spacedim>
ElementValues<spacedim>::~ElementValues()
{
    if (ref_data) delete ref_data;

    for (unsigned int sid=0; sid<n_sides_; sid++)
        for (unsigned int pid=0; pid<n_side_permutations_; pid++)
            delete side_ref_data[sid][pid];
}


template<unsigned int spacedim>
void ElementValues<spacedim>::reinit(const ElementAccessor<spacedim> & cell)
{
	OLD_ASSERT_EQUAL( dim_, cell.dim() );
    data.cell = cell;

    // calculate Jacobian of mapping, JxW, inverse Jacobian
    switch (dim_)
    {
        case 1:
            fill_data<1>();
            break;
        case 2:
            fill_data<2>();
            break;
        case 3:
            fill_data<3>();
            break;
        default:
            ASSERT(false)(dim_).error("Unsupported dimension.\n");
            break;
    }
}


template<unsigned int spacedim>
void ElementValues<spacedim>::reinit(const Side & cell_side)
{
    ASSERT_EQ_DBG( dim_, cell_side.dim() );
    data.side = cell_side;
    
    // calculate Jacobian of mapping, JxW, inverse Jacobian, normal vector(s)
    switch (dim_)
    {
        case 1:
            fill_data<1>();
            fill_side_data<1>();
            break;
        case 2:
            fill_data<2>();
            fill_side_data<2>();
            break;
        case 3:
            fill_data<3>();
            fill_side_data<3>();
            break;
        default:
            ASSERT(false)(dim_).error("Unsupported dimension.\n");
            break;
    }
}


template<unsigned int spacedim>
template<unsigned int dim>
void ElementValues<spacedim>::fill_data()
{
    typename MappingP1<dim,spacedim>::ElementMap coords;

    if ((data.update_flags & update_jacobians) |
        (data.update_flags & update_volume_elements) |
        (data.update_flags & update_JxW_values) |
        (data.update_flags & update_inverse_jacobians) |
        (data.update_flags & update_normal_vectors) |
        (data.update_flags & update_quadrature_points))
    {
        if (cell().is_valid())
            coords = MappingP1<dim,spacedim>::element_map(cell());
        else
            coords = MappingP1<dim,spacedim>::element_map(side().element());
    }

    // calculation of Jacobian dependent data
    if ((data.update_flags & update_jacobians) |
        (data.update_flags & update_volume_elements) |
        (data.update_flags & update_JxW_values) |
        (data.update_flags & update_inverse_jacobians) |
        (data.update_flags & update_normal_vectors))
    {
        arma::mat::fixed<spacedim,dim> jac = MappingP1<dim,spacedim>::jacobian(coords);

        // update Jacobians
        if (data.update_flags & update_jacobians)
            for (unsigned int i=0; i<n_points_; i++)
                data.jacobians.set(i) = Armor::mat<spacedim,dim>( jac );

        // calculation of determinant dependent data
        if (data.update_flags & update_volume_elements | update_JxW_values)
        {
            double det = fabs(::determinant(jac));

            // update determinants
            if (data.update_flags & update_volume_elements)
                for (unsigned int i=0; i<n_points_; i++)
                    data.determinants[i] = det;

            // update JxW values
            if (data.update_flags & update_JxW_values)
                for (unsigned int i=0; i<n_points_; i++)
                    data.JxW_values[i] = det*ref_data->weights[i];
        }

        // update inverse Jacobians
        if (data.update_flags & update_inverse_jacobians)
        {
            arma::mat::fixed<dim,spacedim> ijac;
            if (dim==spacedim)
            {
                ijac = inv(jac);
            }
            else
            {
                ijac = pinv(jac);
            }
            for (unsigned int i=0; i<n_points_; i++)
                data.inverse_jacobians.set(i) = Armor::mat<dim,spacedim>( ijac );
        }
    }

    // quadrature points in the actual cell coordinate system
    if (cell().is_valid() && data.update_flags & update_quadrature_points)
    {
        for (unsigned int i=0; i<n_points_; i++)
            data.points.set(i) = Armor::vec<spacedim>( coords*ref_data->bar_coords[i] );
    }
}


template<unsigned int spacedim>
template<unsigned int dim>
void ElementValues<spacedim>::fill_side_data()
{
    const unsigned int side_idx = side().side_idx();
    const unsigned int perm_idx = side().element()->permutation_idx(side_idx);

    // calculation of normal vectors to the side
    if (data.update_flags & update_normal_vectors)
    {
        arma::vec::fixed<spacedim> n_cell;
        n_cell = trans(data.inverse_jacobians.template mat<dim,spacedim>(0))*RefElement<dim>::normal_vector(side_idx);
        n_cell = n_cell/norm(n_cell,2);
        for (unsigned int i=0; i<n_points_; i++)
            data.normal_vectors.set(i) = Armor::vec<spacedim>( n_cell );
    }

    // Quadrature points in the actual cell coordinate system.
    // The points location can vary from side to side.
    if (data.update_flags & update_quadrature_points)
    {
        typename MappingP1<dim,spacedim>::ElementMap coords = MappingP1<dim,spacedim>::element_map(side().element());
        for (unsigned int i=0; i<n_points_; i++)
            data.points.set(i) = Armor::vec<spacedim>( coords*side_ref_data[side_idx][perm_idx]->bar_coords[i] );
    }

    if (data.update_flags & update_side_JxW_values)
    {
        double side_det;
        if (dim <= 1)
        {
            side_det = 1;
        }
        else
        {
            arma::mat::fixed<spacedim,dim> side_coords;
            arma::mat::fixed<spacedim, MatrixSizes<dim>::dim_minus_one > side_jac;   // some compilers complain for case dim==0

            // calculation of side Jacobian
            for (unsigned int n=0; n<dim; n++)
                for (unsigned int c=0; c<spacedim; c++)
                    side_coords(c,n) = (*data.side.node(n))[c];
            side_jac = MappingP1<MatrixSizes<dim>::dim_minus_one,spacedim>::jacobian(side_coords);

            // calculation of JxW
            side_det = fabs(::determinant(side_jac));
        }
        for (unsigned int i=0; i<n_points_; i++)
            data.side_JxW_values[i] = side_det*side_ref_data[side_idx][perm_idx]->weights[i];
    }
}





template class ElementValues<3>;

