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


template<unsigned int dim, unsigned int spacedim>
void ElementData<dim,spacedim>::allocate(unsigned int size, UpdateFlags flags)
{
    update_flags = flags;

    // resize the arrays of computed quantities
    if (update_flags & update_jacobians)
        jacobians.resize(size);

    if (update_flags & update_volume_elements)
        determinants.resize(size);

    if ((update_flags & update_JxW_values) |
        (update_flags & update_side_JxW_values))
        JxW_values.resize(size);

    if (update_flags & update_inverse_jacobians)
        inverse_jacobians.resize(size);

    if (update_flags & update_quadrature_points)
        points.resize(size);

    if (update_flags & update_normal_vectors)
        normal_vectors.resize(size);
}


template<unsigned int dim, unsigned int spacedim>
void ElementData<dim,spacedim>::print()
{
    if (present_cell.is_valid())
    {
        printf("cell %d dim %d ", present_cell.idx(), present_cell.dim());
        
        printf(" det[");
        for (auto d : determinants) printf("%g ", d); printf("]");
        
        printf(" JxW[");
        for (auto j : JxW_values) printf("%g ", j); printf("]");
        
        printf(" nv[");
        for (auto n : normal_vectors)
        {
            printf(" [");
            for (unsigned int c=0; c<spacedim; c++) printf("%g ", n[c]);
            printf("]");
        }
        printf("]\n");
    }
}


template<unsigned int dim,unsigned int spacedim>
ElementValuesBase<dim,spacedim>::ElementValuesBase()
: n_points_(0)
{}



template<unsigned int dim,unsigned int spacedim>
ElementValuesBase<dim,spacedim>::~ElementValuesBase()
{}



template<unsigned int dim, unsigned int spacedim>
void ElementValuesBase<dim,spacedim>::allocate(unsigned int n_points, UpdateFlags _flags)
{
    n_points_ = n_points;
    data.allocate(n_points_, update_each(_flags));
}



template<unsigned int dim, unsigned int spacedim>
RefElementData *ElementValuesBase<dim,spacedim>::init_ref_data(const Quadrature &q)
{
    ASSERT_DBG( q.dim() == dim );
    RefElementData *data = new RefElementData(q.size());

    for (unsigned int i=0; i<q.size(); i++)
    {
        data->bar_coords[i] = RefElement<dim>::local_to_bary(q.point<dim>(i));
        data->weights[i] = q.weight(i);
    }
    
    return data;
}



template<unsigned int dim, unsigned int spacedim>
UpdateFlags ElementValuesBase<dim,spacedim>::update_each(UpdateFlags flags)
{
    return MappingP1<dim,spacedim>::update_each(flags);
}



template<unsigned int dim, unsigned int spacedim>
ElementValues<dim,spacedim>::ElementValues(
         Quadrature &_quadrature,
         UpdateFlags _flags)
: ElementValuesBase<dim, spacedim>(),
  quadrature_(&_quadrature),
  ref_data(nullptr)
{
    if (dim == 0) return; // avoid unnecessary allocation of dummy 0 dimensional objects
    ASSERT_DBG( _quadrature.dim() == dim );
    this->allocate(_quadrature.size(), _flags);

    // precompute finite element data
    ref_data = this->init_ref_data(_quadrature);
}


template<unsigned int dim, unsigned int spacedim>
ElementValues<dim,spacedim>::~ElementValues()
{
    if (ref_data) delete ref_data;
}



template<unsigned int dim,unsigned int spacedim>
void ElementValues<dim,spacedim>::reinit(const ElementAccessor<3> & cell)
{
	OLD_ASSERT_EQUAL( dim, cell->dim() );
    this->data.present_cell = cell;

    // calculate Jacobian of mapping, JxW, inverse Jacobian
    fill_data();
}


template<unsigned int dim, unsigned int spacedim>
void ElementValues<dim,spacedim>::fill_data()
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
                this->data.jacobians[i] = jac;

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
                this->data.inverse_jacobians[i] = ijac;
        }
    }

    // quadrature points in the actual cell coordinate system
    if (this->data.update_flags & update_quadrature_points)
    {
        for (unsigned int i=0; i<this->n_points_; i++)
            this->data.points[i] = coords*ref_data->bar_coords[i];
    }
}









template<unsigned int dim,unsigned int spacedim>
ElemSideValues<dim,spacedim>::ElemSideValues(
                                 Quadrature & _sub_quadrature,
                                 const UpdateFlags _flags)
: ElementValuesBase<dim,spacedim>(),
  side_quad(RefElement<dim>::n_sides, std::vector<Quadrature>(RefElement<dim>::n_side_permutations, Quadrature(dim)))
{
    ASSERT_DBG( _sub_quadrature.dim() + 1 == dim );
    
    this->allocate(_sub_quadrature.size(), _flags);

    for (unsigned int sid = 0; sid < RefElement<dim>::n_sides; sid++)
    {
    	for (unsigned int pid = 0; pid < RefElement<dim>::n_side_permutations; pid++)
    	{
    		// transform the side quadrature points to the cell quadrature points
            side_quad[sid][pid] = _sub_quadrature.make_from_side<dim>(sid, pid);
    		side_ref_data[sid][pid] = this->init_ref_data(side_quad[sid][pid]);
    	}
    }
}



template<unsigned int dim,unsigned int spacedim>
ElemSideValues<dim,spacedim>::~ElemSideValues()
{
    for (unsigned int sid=0; sid<RefElement<dim>::n_sides; sid++)
        for (unsigned int pid=0; pid<RefElement<dim>::n_side_permutations; pid++)
            delete side_ref_data[sid][pid];
}


template<unsigned int dim,unsigned int spacedim>
void ElemSideValues<dim,spacedim>::reinit(const ElementAccessor<3> & cell,
		unsigned int sid)
{
    ASSERT_LT_DBG( sid, cell->n_sides() );
    ASSERT_EQ_DBG( dim, cell->dim() );
    this->data.present_cell = cell;
    side_idx_ = sid;
    
    // calculate Jacobian of mapping, JxW, inverse Jacobian, normal vector(s)
    fill_data();
}


template<unsigned int dim, unsigned int spacedim>
void ElemSideValues<dim,spacedim>::fill_data()
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
                this->data.jacobians[i] = jac;

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
            ASSERT_LE_DBG(this->n_points_, this->data.inverse_jacobians.size());
            for (unsigned int i=0; i<this->n_points_; i++)
                this->data.inverse_jacobians[i] = ijac;

            // calculation of normal vectors to the side
            if ((this->data.update_flags & update_normal_vectors))
            {
                arma::vec::fixed<spacedim> n_cell;
                n_cell = trans(ijac)*RefElement<dim>::normal_vector(side_idx_);
                n_cell = n_cell/norm(n_cell,2);
                for (unsigned int i=0; i<this->n_points_; i++)
                    this->data.normal_vectors[i] = n_cell;
            }
        }
    }

    // Quadrature points in the actual cell coordinate system.
    // The points location can vary from side to side.
    if (this->data.update_flags & update_quadrature_points)
    {
        for (unsigned int i=0; i<this->n_points_; i++)
            this->data.points[i] = coords*side_ref_data[side_idx_][this->data.present_cell->permutation_idx(side_idx_)]->bar_coords[i];
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
                    side_coords(c,n) = (*this->data.present_cell.side(side_idx_)->node(n))[c];
            side_jac = MappingP1<MatrixSizes<dim>::dim_minus_one,spacedim>::jacobian(side_coords);

            // calculation of JxW
            side_det = fabs(determinant(side_jac));
        }
        for (unsigned int i=0; i<this->n_points_; i++)
            this->data.JxW_values[i] = side_det*side_ref_data[side_idx_][this->data.present_cell->permutation_idx(side_idx_)]->weights[i];
    }
}





template class ElementData<0,3>;
template class ElementData<1,3>;
template class ElementData<2,3>;
template class ElementData<3,3>;

template class ElementValuesBase<0,3>;
template class ElementValuesBase<1,3>;
template class ElementValuesBase<2,3>;
template class ElementValuesBase<3,3>;

template class ElementValues<0,3>;
template class ElementValues<1,3>;
template class ElementValues<2,3>;
template class ElementValues<3,3>;

template class ElemSideValues<1,3>;
template class ElemSideValues<2,3>;
template class ElemSideValues<3,3>;

