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
 * @brief Abstract class for description of finite elements.
 * @author Jan Stebel
 */


#include "system/system.hh"
#include "quadrature/quadrature.hh"
#include "fem/dofhandler.hh"
#include "fem/finite_element.hh"
#include "fem/fe_values.hh"



using namespace std;




template<unsigned int dim, unsigned int spacedim>
FiniteElement<dim,spacedim>::FiniteElement()
{
    init();
}

template<unsigned int dim, unsigned int spacedim>
void FiniteElement<dim,spacedim>::init()
{
    number_of_dofs = 0;
    is_scalar_fe = true;
    for (int i = 0; i <= dim; i++)
    {
        number_of_single_dofs[i] = 0;
        number_of_pairs[i] = 0;
        number_of_triples[i] = 0;
        number_of_sextuples[i] = 0;
    }
}

template<unsigned int dim, unsigned int spacedim> inline
const unsigned int FiniteElement<dim,spacedim>::n_dofs()
{
    return number_of_dofs;
}

template<unsigned int dim, unsigned int spacedim> inline
const unsigned int FiniteElement<dim,spacedim>::n_object_dofs(
        unsigned int object_dim, DofMultiplicity multiplicity)
{
    ASSERT(object_dim >= 0 && object_dim <= dim,
            "Object type number is out of range.");
    switch (multiplicity)
    {
    case DOF_SINGLE:
        return number_of_single_dofs[object_dim];
    case DOF_PAIR:
        return number_of_pairs[object_dim];
    case DOF_TRIPLE:
        return number_of_triples[object_dim];
    case DOF_SEXTUPLE:
        return number_of_sextuples[object_dim];
    }
}

template<unsigned int dim, unsigned int spacedim> inline
void FiniteElement<dim,spacedim>::compute_node_matrix()
{
    ASSERT_SIZES(get_generalized_support_points().size(), number_of_dofs);

    arma::mat M(number_of_dofs, number_of_dofs);

    for (int i = 0; i < number_of_dofs; i++)
        for (int j = 0; j < number_of_dofs; j++) {
            M(j, i) = basis_value(j, get_generalized_support_points()[i]);

        }
    node_matrix = arma::inv(M);
}

template<unsigned int dim, unsigned int spacedim>
FEInternalData *FiniteElement<dim,spacedim>::initialize(const Quadrature<dim> &q, UpdateFlags flags)
{
    FEInternalData *data = new FEInternalData;

    if (flags & update_values)
    {
        arma::vec values(number_of_dofs);
        data->basis_values.resize(q.size());
        for (int i=0; i<q.size(); i++)
        {
            for (int j=0; j<number_of_dofs; j++)
                values[j] = basis_value(j, q.point(i));
            data->basis_values[i] = node_matrix * values;
        }
    }

    if (flags & update_gradients)
    {
        arma::mat grads(number_of_dofs, dim);
        data->basis_grads.resize(q.size());
        for (int i=0; i<q.size(); i++)
        {
            for (int j=0; j<number_of_dofs; j++)
                grads.row(j) = arma::trans(basis_grad(j, q.point(i)));
            data->basis_grads[i] = node_matrix * grads;
        }
    }

    return data;
}

template<unsigned int dim, unsigned int spacedim> inline
UpdateFlags FiniteElement<dim,spacedim>::update_each(UpdateFlags flags)
{
    UpdateFlags f = flags;

    if (flags & update_gradients)
        f |= update_inverse_jacobians;

    return f;
}

template<unsigned int dim, unsigned int spacedim> inline
void FiniteElement<dim,spacedim>::fill_fe_values(
        const Quadrature<dim> &q,
        FEInternalData &data,
        FEValuesData<dim,spacedim> &fv_data)
{
    // shape values
    if (fv_data.update_flags & update_values)
    {
        for (int i = 0; i < q.size(); i++)
            fv_data.shape_values[i] = data.basis_values[i];
    }

    // shape gradients
    if (fv_data.update_flags & update_gradients)
    {
        for (int i = 0; i < q.size(); i++)
        {
            fv_data.shape_gradients[i] = data.basis_grads[i] * fv_data.inverse_jacobians[i];
        }
    }
}

template<unsigned int dim, unsigned int spacedim>
const vector<arma::vec::fixed<dim> > &FiniteElement<dim,spacedim>::get_generalized_support_points()
{
    if (generalized_support_points.size() > 0)
    {
        return generalized_support_points;
    }
    else
    {
        return unit_support_points;
    }
}


template class FiniteElement<0,3>;
template class FiniteElement<1,3>;
template class FiniteElement<2,3>;
template class FiniteElement<3,3>;


