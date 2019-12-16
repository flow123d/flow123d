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
#include "mesh/side_impl.hh"



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
  JxW_values(flags & update_JxW_values | update_side_JxW_values ? size : 0),
  jacobians(flags & update_jacobians ? size : 0, spacedim, dim),
  determinants(flags & update_volume_elements ? size : 0),
  inverse_jacobians(flags & update_inverse_jacobians ? size : 0, dim, spacedim),
  points(flags & update_quadrature_points ? size : 0, spacedim),
  normal_vectors(flags & update_normal_vectors ? size : 0, spacedim),
  update_flags(flags)
{}


template<unsigned int spacedim>
void ElementData<spacedim>::print()
{
    if (present_cell.is_valid())
    {
        printf("cell %d dim %d ", present_cell.idx(), present_cell.dim());
        
        printf(" det[");
        for (auto d : determinants) printf("%g ", d); printf("]");
        
        printf(" JxW[");
        for (auto j : JxW_values) printf("%g ", j); printf("]");
        
        printf(" nv[");
        for (unsigned int i=0; i<normal_vectors.n_vals(); i++)
        {
            auto n = normal_vectors.arma_vec(i);
            printf(" [");
            for (unsigned int c=0; c<spacedim; c++) printf("%g ", n[c]);
            printf("]");
        }
        printf("]\n");
    }
}





template<unsigned int spacedim>
RefElementData *ElementValuesBase<spacedim>::init_ref_data(const Quadrature &q)
{
    ASSERT_DBG( q.dim() == this->data.dim_ );
    ASSERT_DBG( q.size() == n_points_ );
    RefElementData *ref_data = new RefElementData(q.size());

    for (unsigned int i=0; i<q.size(); i++)
    {
        switch (q.dim())
        {
            case 1:
                ref_data->bar_coords[i] = RefElement<1>::local_to_bary(q.point<1>(i).arma());
                break;
            case 2:
                ref_data->bar_coords[i] = RefElement<2>::local_to_bary(q.point<2>(i).arma());
                break;
            case 3:
                ref_data->bar_coords[i] = RefElement<3>::local_to_bary(q.point<3>(i).arma());
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
UpdateFlags ElementValuesBase<spacedim>::update_each(UpdateFlags flags)
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
: ElementValuesBase<spacedim>(_quadrature.size(), _flags, dim ),
  quadrature_(&_quadrature),
  ref_data(nullptr)
{
    if (dim == 0) return; // avoid unnecessary allocation of dummy 0 dimensional objects
    ASSERT_DBG( _quadrature.dim() == dim );

    // precompute finite element data
    ref_data = this->init_ref_data(_quadrature);
}


template<unsigned int spacedim>
ElementValues<spacedim>::~ElementValues()
{
    if (ref_data) delete ref_data;
}



template<unsigned int spacedim>
void ElementValues<spacedim>::reinit(const ElementAccessor<spacedim> & cell)
{
	OLD_ASSERT_EQUAL( this->dim_, cell->dim() );
    this->data.present_cell = cell;

    // calculate Jacobian of mapping, JxW, inverse Jacobian
    switch (this->dim_)
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
            ASSERT(false)(this->dim_).error("Unsupported dimension.\n");
            break;
    }
}


template<unsigned int spacedim>
template<unsigned int dim>
void ElementValues<spacedim>::fill_data()
{
    typename MappingP1<dim,spacedim>::ElementMap coords;
    arma::mat::fixed<spacedim,dim> jac;

    if ((this->data.update_flags & update_jacobians) |
        (this->data.update_flags & update_volume_elements) |
        (this->data.update_flags & update_JxW_values) |
        (this->data.update_flags & update_inverse_jacobians) |
        (this->data.update_flags & update_quadrature_points))
    {
        coords = MappingP1<dim,spacedim>::element_map(this->data.present_cell);
    }

    // calculation of Jacobian dependent data
    if ((this->data.update_flags & update_jacobians) |
        (this->data.update_flags & update_volume_elements) |
        (this->data.update_flags & update_JxW_values) |
        (this->data.update_flags & update_inverse_jacobians))
    {
        jac = MappingP1<dim,spacedim>::jacobian(coords);

        // update Jacobians
        if (this->data.update_flags & update_jacobians)
            for (unsigned int i=0; i<this->n_points_; i++)
                this->data.jacobians.get<spacedim,dim>(i) = jac;

        // calculation of determinant dependent data
        if ((this->data.update_flags & update_volume_elements) | (this->data.update_flags & update_JxW_values))
        {
            double det = fabs(determinant(jac));

            // update determinants
            if (this->data.update_flags & update_volume_elements)
                for (unsigned int i=0; i<this->n_points_; i++)
                    this->data.determinants[i] = det;

            // update JxW values
            if (this->data.update_flags & update_JxW_values)
                for (unsigned int i=0; i<this->n_points_; i++)
                    this->data.JxW_values[i] = det*ref_data->weights[i];
        }

        // update inverse Jacobians
        if (this->data.update_flags & update_inverse_jacobians)
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
            for (unsigned int i=0; i<this->n_points_; i++)
                this->data.inverse_jacobians.get<dim,spacedim>(i) = ijac;
        }
    }

    // quadrature points in the actual cell coordinate system
    if (this->data.update_flags & update_quadrature_points)
    {
        for (unsigned int i=0; i<this->n_points_; i++)
            this->data.points.get<spacedim>(i) = coords*ref_data->bar_coords[i];
    }
}









template<unsigned int spacedim>
ElemSideValues<spacedim>::ElemSideValues(
                                 Quadrature & _sub_quadrature,
                                 const UpdateFlags _flags,
                                 unsigned int dim)
: ElementValuesBase<spacedim>(_sub_quadrature.size(), _flags, dim),
  n_sides_(dim+1),
  n_side_permutations_((dim+1)*(2*dim*dim-5*dim+6)/6),
  side_quad(n_sides_, std::vector<Quadrature>(n_side_permutations_, Quadrature(dim))),
  side_ref_data(n_sides_, std::vector<RefElementData*>(n_side_permutations_))
{
    ASSERT_DBG( _sub_quadrature.dim() + 1 == dim );

    for (unsigned int sid = 0; sid < n_sides_; sid++)
    {
    	for (unsigned int pid = 0; pid < n_side_permutations_; pid++)
    	{
    		// transform the side quadrature points to the cell quadrature points
            switch (dim)
            {
                case 1:
                    side_quad[sid][pid] = _sub_quadrature.make_from_side<1>(sid, pid);
                    break;
                case 2:
                    side_quad[sid][pid] = _sub_quadrature.make_from_side<2>(sid, pid);
                    break;
                case 3:
                    side_quad[sid][pid] = _sub_quadrature.make_from_side<3>(sid, pid);
                    break;
                default:
                    ASSERT(false)(dim).error("Unsupported dimension.\n");
                    break;
            }
    		side_ref_data[sid][pid] = this->init_ref_data(side_quad[sid][pid]);
    	}
    }
}



template<unsigned int spacedim>
ElemSideValues<spacedim>::~ElemSideValues()
{
    for (unsigned int sid=0; sid<n_sides_; sid++)
        for (unsigned int pid=0; pid<n_side_permutations_; pid++)
            delete side_ref_data[sid][pid];
}


template<unsigned int spacedim>
void ElemSideValues<spacedim>::reinit(const ElementAccessor<spacedim> & cell,
		unsigned int sid)
{
    ASSERT_LT_DBG( sid, cell->n_sides() );
    ASSERT_EQ_DBG( this->dim_, cell->dim() );
    this->data.present_cell = cell;
    side_idx_ = sid;
    
    // calculate Jacobian of mapping, JxW, inverse Jacobian, normal vector(s)
    switch (this->dim_)
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
            ASSERT(false)(this->dim_).error("Unsupported dimension.\n");
            break;
    }
}


template<unsigned int spacedim>
template<unsigned int dim>
void ElemSideValues<spacedim>::fill_data()
{
    typename MappingP1<dim,spacedim>::ElementMap coords;

    if ((this->data.update_flags & update_jacobians) |
        (this->data.update_flags & update_volume_elements) |
        (this->data.update_flags & update_inverse_jacobians) |
        (this->data.update_flags & update_normal_vectors) |
        (this->data.update_flags & update_quadrature_points))
    {
        coords = MappingP1<dim,spacedim>::element_map(this->data.present_cell);
    }

    // calculation of cell Jacobians and dependent data
    if ((this->data.update_flags & update_jacobians) |
        (this->data.update_flags & update_volume_elements) |
        (this->data.update_flags & update_inverse_jacobians) |
        (this->data.update_flags & update_normal_vectors))
    {
        arma::mat::fixed<spacedim,dim> jac = MappingP1<dim,spacedim>::jacobian(coords);

        // update cell Jacobians
        if (this->data.update_flags & update_jacobians)
            for (unsigned int i=0; i<this->n_points_; i++)
                this->data.jacobians.get<spacedim,dim>(i) = jac;

        // update determinants of Jacobians
        if (this->data.update_flags & update_volume_elements)
        {
            double det = fabs(determinant(jac));
            for (unsigned int i=0; i<this->n_points_; i++)
                this->data.determinants[i] = det;
        }

        // inverse Jacobians
        if (this->data.update_flags & update_inverse_jacobians)
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
            ASSERT_LE_DBG(this->n_points_, this->data.inverse_jacobians.n_vals());
            for (unsigned int i=0; i<this->n_points_; i++)
                this->data.inverse_jacobians.get<dim,spacedim>(i) = ijac;

            // calculation of normal vectors to the side
            if ((this->data.update_flags & update_normal_vectors))
            {
                arma::vec::fixed<spacedim> n_cell;
                n_cell = trans(ijac)*RefElement<dim>::normal_vector(side_idx_);
                n_cell = n_cell/norm(n_cell,2);
                for (unsigned int i=0; i<this->n_points_; i++)
                    this->data.normal_vectors.get<spacedim>(i) = n_cell;
            }
        }
    }

    // Quadrature points in the actual cell coordinate system.
    // The points location can vary from side to side.
    if (this->data.update_flags & update_quadrature_points)
    {
        for (unsigned int i=0; i<this->n_points_; i++)
            this->data.points.get<spacedim>(i) = coords*side_ref_data[side_idx_][this->data.present_cell->permutation_idx(side_idx_)]->bar_coords[i];
    }

    if (this->data.update_flags & update_side_JxW_values)
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
                    side_coords(c,n) = this->data.present_cell.side(side_idx_)->node(n)->point()[c];
            side_jac = MappingP1<MatrixSizes<dim>::dim_minus_one,spacedim>::jacobian(side_coords);

            // calculation of JxW
            side_det = fabs(determinant(side_jac));
        }
        for (unsigned int i=0; i<this->n_points_; i++)
            this->data.JxW_values[i] = side_det*side_ref_data[side_idx_][this->data.present_cell->permutation_idx(side_idx_)]->weights[i];
    }
}





template class ElementValuesBase<3>;
template class ElementValues<3>;
template class ElemSideValues<3>;

