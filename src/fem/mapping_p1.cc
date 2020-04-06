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
 * @file    mapping_p1.cc
 * @brief   Class MappingP1 implements the affine transformation of
 *          the unit cell onto the actual cell.
 * @author  Jan Stebel
 */

#include "fem/mapping_p1.hh"
#include "quadrature/quadrature.hh"
#include "fem/fe_values.hh"



using namespace std;



template<unsigned int dim, unsigned int spacedim>
MappingP1<dim,spacedim>::MappingP1()
{
}

template<unsigned int dim, unsigned int spacedim>
MappingInternalData *MappingP1<dim,spacedim>::initialize(const Quadrature &q, UpdateFlags flags)
{
    ASSERT_DBG( q.dim() == dim );
    MappingInternalData *data = new MappingInternalData;

    // Initialize the gradients of the canonical basis in the
    // barycentric coordinates.
    // In the case of P1 mapping the shape functions are linear,
    // hence the gradients are constant and can be precomputed.
    if ((flags & update_jacobians) |
            (flags & update_volume_elements) |
            (flags & update_JxW_values) |
            (flags & update_side_JxW_values) |
            (flags & update_inverse_jacobians))
    {
        grad.zeros();
        for (unsigned int i=0; i<dim; i++)
        {
            grad(0,i) = -1;
            grad(i+1,i) = 1;
        }
    }

    // barycentric coordinates of quadrature points
    if (flags & update_quadrature_points)
    {
        data->bar_coords.resize(q.size());
        for (unsigned int i=0; i<q.size(); i++)
            data->bar_coords[i] = RefElement<dim>::local_to_bary(q.point<dim>(i).arma());
    }



    return data;
}

template<unsigned int dim, unsigned int spacedim>
UpdateFlags MappingP1<dim,spacedim>::update_each(UpdateFlags flags)
{
    UpdateFlags f = flags;
    
    if (flags & update_normal_vectors)
        f |= update_inverse_jacobians;

    if ((flags & update_volume_elements) |
        (flags & update_JxW_values) |
        (flags & update_inverse_jacobians))
        f |= update_jacobians;

    return f;
}


template<unsigned int dim, unsigned int spacedim>
void MappingP1<dim,spacedim>::fill_fe_values(const ElementAccessor<3> &cell,
                            const Quadrature &q,
                            MappingInternalData &data,
                            FEValuesData<dim,spacedim> &fv_data)
{
    ASSERT_DBG( q.dim() == dim );
    ElementMap coords;
    arma::mat::fixed<spacedim,dim> jac;

    if ((fv_data.update_flags & update_jacobians) |
        (fv_data.update_flags & update_volume_elements) |
        (fv_data.update_flags & update_JxW_values) |
        (fv_data.update_flags & update_inverse_jacobians) |
        (fv_data.update_flags & update_quadrature_points))
    {
        coords = element_map(cell);
    }

    // calculation of Jacobian dependent data
    if ((fv_data.update_flags & update_jacobians) |
        (fv_data.update_flags & update_volume_elements) |
        (fv_data.update_flags & update_JxW_values) |
        (fv_data.update_flags & update_inverse_jacobians))
    {
        jac = coords*grad;

        // update Jacobians
        if (fv_data.update_flags & update_jacobians)
            for (unsigned int i=0; i<q.size(); i++)
                fv_data.jacobians[i] = jac;

        // calculation of determinant dependent data
        if ((fv_data.update_flags & update_volume_elements) | (fv_data.update_flags & update_JxW_values))
        {
            double det = fabs(determinant(jac));

            // update determinants
            if (fv_data.update_flags & update_volume_elements)
                for (unsigned int i=0; i<q.size(); i++)
                    fv_data.determinants[i] = det;

            // update JxW values
            if (fv_data.update_flags & update_JxW_values)
                for (unsigned int i=0; i<q.size(); i++)
                    fv_data.JxW_values[i] = det*q.weight(i);
        }

        // update inverse Jacobians
        if (fv_data.update_flags & update_inverse_jacobians)
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
            for (unsigned int i=0; i<q.size(); i++)
                fv_data.inverse_jacobians[i] = ijac;
        }
    }

    // quadrature points in the actual cell coordinate system
    if (fv_data.update_flags & update_quadrature_points)
    {
        BaryPoint basis;
        for (unsigned int i=0; i<q.size(); i++)
            fv_data.points[i] = coords*data.bar_coords[i];
    }
}

template<unsigned int dim, unsigned int spacedim>
void MappingP1<dim,spacedim>::fill_fe_side_values(const ElementAccessor<3> &cell,
                            unsigned int sid,
                            const Quadrature &q,
                            MappingInternalData &data,
                            FEValuesData<dim,spacedim> &fv_data)
{
    ASSERT_DBG( q.dim() == dim );
    ElementMap coords;

    if ((fv_data.update_flags & update_jacobians) |
        (fv_data.update_flags & update_volume_elements) |
        (fv_data.update_flags & update_inverse_jacobians) |
        (fv_data.update_flags & update_normal_vectors) |
        (fv_data.update_flags & update_quadrature_points))
    {
        coords = element_map(cell);
    }

    // calculation of cell Jacobians and dependent data
    if ((fv_data.update_flags & update_jacobians) |
        (fv_data.update_flags & update_volume_elements) |
        (fv_data.update_flags & update_inverse_jacobians) |
        (fv_data.update_flags & update_normal_vectors))
    {
        arma::mat::fixed<spacedim,dim> jac = coords*grad;

        // update cell Jacobians
        if (fv_data.update_flags & update_jacobians)
            for (unsigned int i=0; i<q.size(); i++)
                fv_data.jacobians[i] = jac;

        // update determinants of Jacobians
        if (fv_data.update_flags & update_volume_elements)
        {
            double det = fabs(determinant(jac));
            for (unsigned int i=0; i<q.size(); i++)
                fv_data.determinants[i] = det;
        }

        // inverse Jacobians
        if (fv_data.update_flags & update_inverse_jacobians)
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
            ASSERT_LE_DBG(q.size(), fv_data.inverse_jacobians.size());
            for (unsigned int i=0; i<q.size(); i++)
                fv_data.inverse_jacobians[i] = ijac;

            // calculation of normal vectors to the side
            if ((fv_data.update_flags & update_normal_vectors))
            {
                arma::vec::fixed<spacedim> n_cell;
                n_cell = trans(ijac)*RefElement<dim>::normal_vector(sid);
                n_cell = n_cell/norm(n_cell,2);
                for (unsigned int i=0; i<q.size(); i++)
                    fv_data.normal_vectors[i] = n_cell;
            }
        }
    }

    // Quadrature points in the actual cell coordinate system.
    // The points location can vary from side to side.
    if (fv_data.update_flags & update_quadrature_points)
    {
        for (unsigned int i=0; i<q.size(); i++)
            fv_data.points[i] = coords*data.bar_coords[i];
    }

    if (fv_data.update_flags & update_side_JxW_values)
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
            side_coords.zeros();
            for (unsigned int n=0; n<dim; n++)
                for (unsigned int c=0; c<spacedim; c++)
                    side_coords(c,n) = cell.side(sid)->node(n)->point()[c];
            side_jac = side_coords * grad.submat(0,0,dim-1,dim-2);

            // calculation of JxW
            side_det = fabs(determinant(side_jac));
        }
        for (unsigned int i=0; i<q.size(); i++)
            fv_data.JxW_values[i] = side_det*q.weight(i);
    }
}


template<unsigned int dim, unsigned int spacedim>
auto MappingP1<dim,spacedim>::element_map(ElementAccessor<3> elm) const -> ElementMap
{
    ElementMap coords;
    for (unsigned int i=0; i<dim+1; i++)
        coords.col(i) = elm.node(i)->point();
    return coords;
}


template<unsigned int dim, unsigned int spacedim>
auto MappingP1<dim,spacedim>::project_real_to_unit(const RealPoint &point, const ElementMap &map) const -> BaryPoint
{
    arma::mat::fixed<3, dim> A = map.cols(1,dim);
    for(unsigned int i=0; i < dim; i++ ) {
        A.col(i) -= map.col(0);
    }
    
    arma::mat::fixed<dim, dim> AtA = A.t()*A;
    arma::vec::fixed<dim> Atb = A.t()*(point - map.col(0));
    arma::vec::fixed<dim+1> bary_coord;
    bary_coord.rows(1, dim) = arma::solve(AtA, Atb);
    bary_coord( 0 ) = 1.0 - arma::sum( bary_coord.rows(1,dim) );
    return bary_coord;
}

template<unsigned int dim, unsigned int spacedim>
auto MappingP1<dim,spacedim>::project_unit_to_real(const BaryPoint &point, const ElementMap &map) const -> RealPoint
{
    return map * point;
}

template<unsigned int dim, unsigned int spacedim>
auto MappingP1<dim,spacedim>::clip_to_element(BaryPoint &barycentric) -> BaryPoint{
    return RefElement<dim>::clip(barycentric);
}

template <unsigned int dim, unsigned int spacedim>
bool MappingP1<dim,spacedim>::contains_point(arma::vec point, ElementAccessor<3> elm)
{
	arma::vec projection = this->project_real_to_unit(point, this->element_map(elm));
	return (projection.min() >= -BoundingBox::epsilon);
}



template class MappingP1<1,3>;
template class MappingP1<2,3>;
template class MappingP1<3,3>;




