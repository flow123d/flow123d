/*!
 *
 * Copyright (C) 2007 Technical University of Liberec.  All rights reserved.
 *
 * Please make a following refer to Flow123d on your project site if you use the program for any purpose,
 * especially for academic research:
 * Flow123d, Research Centre: Advanced Remedial Technologies, Technical University of Liberec, Czech Republic
 *
 * This program is free software; you can redistribute it and/or modify it under the terms
 * of the GNU General Public License version 3 as published by the Free Software Foundation.
 *
 * This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
 * without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 * See the GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License along with this program; if not,
 * write to the Free Software Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA 021110-1307, USA.
 *
 *
 * $Id$
 * $Revision$
 * $LastChangedBy$
 * $LastChangedDate$
 *
 * @file
 * @brief Class MappingP1 implements the affine transformation of
 *        the unit cell onto the actual cell.
 * @author Jan Stebel
 */


#include "fem/mapping_p1.hh"
#include "quadrature/quadrature.hh"
#include "fem/fe_values.hh"



using namespace std;
using namespace arma;




template<unsigned int dim, unsigned int spacedim>
MappingP1<dim,spacedim>::MappingP1()
{
}

template<unsigned int dim, unsigned int spacedim>
MappingInternalData *MappingP1<dim,spacedim>::initialize(const Quadrature<dim> &q, UpdateFlags flags)
{
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
        vec::fixed<dim+1> basis;
        data->bar_coords.resize(q.size());
        for (unsigned int i=0; i<q.size(); i++)
        {
            basis[0] = 1;
            for (unsigned int j=0; j<dim; j++)
            {
                basis[0] -= q.point(i)[j];
                basis[j+1] = q.point(i)[j];
            }
            data->bar_coords[i] = basis;
        }
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
void MappingP1<dim,spacedim>::fill_fe_values(const typename DOFHandlerBase::CellIterator &cell,
                            const Quadrature<dim> &q,
                            MappingInternalData &data,
                            FEValuesData<dim,spacedim> &fv_data)
{
    mat::fixed<spacedim,dim+1> coords;
    mat::fixed<spacedim,dim> jac;

    if ((fv_data.update_flags & update_jacobians) |
        (fv_data.update_flags & update_volume_elements) |
        (fv_data.update_flags & update_JxW_values) |
        (fv_data.update_flags & update_inverse_jacobians) |
        (fv_data.update_flags & update_quadrature_points))
    {
        coords.zeros();
        for (unsigned int n=0; n<dim+1; n++)
            for (unsigned int c=0; c<spacedim; c++)
                coords(c,n) = cell->node[n]->point()[c];
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
            mat::fixed<dim,spacedim> ijac;
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
        vec::fixed<dim+1> basis;
        for (unsigned int i=0; i<q.size(); i++)
            fv_data.points[i] = coords*data.bar_coords[i];
    }
}

template<unsigned int dim, unsigned int spacedim>
void MappingP1<dim,spacedim>::fill_fe_side_values(const typename DOFHandlerBase::CellIterator &cell,
                            unsigned int sid,
                            const Quadrature<dim> &q,
                            MappingInternalData &data,
                            FEValuesData<dim,spacedim> &fv_data)
{
    mat::fixed<spacedim,dim+1> coords;

    if ((fv_data.update_flags & update_jacobians) |
        (fv_data.update_flags & update_volume_elements) |
        (fv_data.update_flags & update_inverse_jacobians) |
        (fv_data.update_flags & update_normal_vectors) |
        (fv_data.update_flags & update_quadrature_points))
    {
        coords.zeros();
        for (unsigned int n=0; n<dim+1; n++)
            for (unsigned int c=0; c<spacedim; c++)
                coords(c,n) = cell->node[n]->point()[c];
    }

    // calculation of cell Jacobians and dependent data
    if ((fv_data.update_flags & update_jacobians) |
        (fv_data.update_flags & update_volume_elements) |
        (fv_data.update_flags & update_inverse_jacobians) |
        (fv_data.update_flags & update_normal_vectors))
    {
        mat::fixed<spacedim,dim> jac = coords*grad;

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
            mat::fixed<dim,spacedim> ijac;
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

            // calculation of normal vectors to the side
            if ((fv_data.update_flags & update_normal_vectors))
            {
                vec::fixed<spacedim> n_cell;
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
            mat::fixed<spacedim,dim> side_coords;
            mat::fixed<spacedim, MatrixSizes<dim>::dim_minus_one > side_jac;   // some compilers complain for case dim==0

            // calculation of side Jacobian
            side_coords.zeros();
            for (unsigned int n=0; n<dim; n++)
                for (unsigned int c=0; c<spacedim; c++)
                    side_coords(c,n) = cell->side(sid)->node(n)->point()[c];
            side_jac = side_coords * grad.submat(0,0,dim-1,dim-2);

            // calculation of JxW
            side_det = fabs(determinant(side_jac));
        }
        for (unsigned int i=0; i<q.size(); i++)
            fv_data.JxW_values[i] = side_det*q.weight(i);
    }
}

template<>
void MappingP1<0,3>::fill_fe_side_values(const DOFHandlerBase::CellIterator &cell,
                            unsigned int sid,
                            const Quadrature<0> &q,
                            MappingInternalData &data,
                            FEValuesData<0,3> &fv_data)
{}



template class MappingP1<0,3>;
template class MappingP1<1,3>;
template class MappingP1<2,3>;
template class MappingP1<3,3>;




