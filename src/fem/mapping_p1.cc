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






template<unsigned int dim, unsigned int spacedim>
MappingP1<dim,spacedim>::MappingP1()
{
    // Initialize the gradients of shape functions.
    // In the case of P1 mapping the shape functions are linear,
    // hence the gradients are constant and can be precomputed.
    grad.zeros();
    for (int i=0; i<dim; i++)
    {
        grad(i,0) = -1;
        grad(i,i+1) = 1;
    }
}


template<unsigned int dim, unsigned int spacedim>
void MappingP1<dim,spacedim>::fill_fe_values(const typename DOFHandler<dim>::CellIterator &cell,
                            const Quadrature<dim> &q,
                            vector< mat::fixed<dim,spacedim> > &jacobians,
                            vector<double> &JxW_values,
                            vector< mat > &inverse_jacobians,
                            vector< vec::fixed<spacedim> > &normal_vectors)
{
    mat::fixed<dim+1,spacedim> coords;
    mat::fixed<dim,spacedim> jac;

    coords.zeros();

    // calculation of Jacobians
    for (int n=0; n<dim+1; n++)
    {
        for (int c=0; c<spacedim; c++)
        {
//            cell->node[n]->point().print();
            double x = (cell->node[n]->point())[c];
            coords(n,c) = x;
        }
    }
    jac = grad*coords;
    for (int i=0; i<q.size(); i++)
        jacobians[i] = jac;

    // calculation of JxW
    double det = determinant(jacobians[0]);
    for (int i=0; i<q.size(); i++)
        JxW_values[i] = det*q.weight(i);

    // inverse Jacobians
    if (dim==spacedim)
    {
        inverse_jacobians[0] = inv(jacobians[0]);
    }
    else
    {
        inverse_jacobians[0] = pinv(jacobians[0]);
    }
    for (int i=1; i<q.size(); i++)
        inverse_jacobians[i] = inverse_jacobians[0];


    // normal vectors...
}



template class MappingP1<2,2>;




