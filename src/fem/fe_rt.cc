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
 * @file    fe_rt.cc
 * @brief   Definitions of Raviart-Thomas finite elements.
 * @author  Jan Stebel
 */

#include "fem/fe_rt.hh"
#include "fem/fe_values.hh"
#include "mesh/ref_element.hh"
#include "quadrature/quadrature_lib.hh"


template<unsigned int dim, unsigned int spacedim>
FE_RT0<dim,spacedim>::FE_RT0()
{
    arma::vec::fixed<dim> sp;

    this->init(spacedim, false, FEVector);
    number_of_dofs = dim+1;
    number_of_single_dofs[dim] = dim+1;
    
    this->component_indices_.clear();
    this->nonzero_components_.resize(number_of_dofs, std::vector<bool>(spacedim, true));

    for (unsigned int sid=0; sid<RefElement<dim>::n_sides; ++sid)
    {
        sp.fill(0);
        for (unsigned int i=0; i<RefElement<dim>::n_nodes_per_side; ++i)
            sp += RefElement<dim>::node_coords(RefElement<dim>::interact(Interaction<0,dim-1>(sid))[i]);
        sp /= RefElement<dim>::n_nodes_per_side;
        generalized_support_points.push_back(sp);
    }

    compute_node_matrix();
}


template<unsigned int dim, unsigned int spacedim>
double FE_RT0<dim,spacedim>::basis_value(const unsigned int i,
            const arma::vec::fixed<dim> &p, const unsigned int comp) const
{
  OLD_ASSERT(i<n_raw_functions, "Index of basis function is out of range.");
  OLD_ASSERT(comp<this->n_components_, "Index of component is out of range.");
  
  if (i>0 && comp==i-1)
    return p[comp]-1;
  else
    return p[comp];
}


template<unsigned int dim, unsigned int spacedim>
arma::vec::fixed<dim> FE_RT0<dim,spacedim>::basis_grad(const unsigned int i,
            const arma::vec::fixed<dim> &p, const unsigned int comp) const
{
  OLD_ASSERT(i<n_raw_functions, "Index of basis function is out of range.");
  OLD_ASSERT(comp<this->n_components_, "Index of component is out of range.");
  
  arma::vec::fixed<dim> g;
  g.zeros();
  g[comp] = 1;

  return g;
}

template<unsigned int dim, unsigned int spacedim>
void FE_RT0<dim,spacedim>::compute_node_matrix()
{
    arma::mat::fixed<n_raw_functions,dim+1> F;
    arma::vec::fixed<dim> r;

    /*
     * Node matrix helps creating the shape functions $\{b_k\}$ from
     * the raw basis $\{r_i\}$:
     *
     * $$  b_k = \sum_{i=0}^{dim*(dim+1)-1} N_{ki} r_i. $$
     *
     * The shape functions must obey the flux condition
     *
     * $$ b_k\cdot n_j |\Gamma_j| = \delta_{kj}, $$
     *
     * where $n_j$, $|\Gamma_j|$ is the unit outward normal vector and
     * the area of the $j$-th side, respectively. Consequently,
     * the node matrix $N$ is determined as the Moon-Penrose
     * pseudoinverse of the flux matrix, i.e.:
     *
     * $$ NF = I,\quad F_{ij} = r_i\cdot n_j |\Gamma_j|. $$
     *
     */

    for (unsigned int i=0; i<n_raw_functions; i++)
    {
        for (unsigned int j=0; j<dim+1; ++j)
        {
            double dot_product = 0;
            for (unsigned int c=0; c<dim; c++)
               dot_product += basis_value(i,generalized_support_points[j],c)*RefElement<dim>::normal_vector(j)(c);
            F(i,j) = dot_product*RefElement<dim>::side_measure(j);
        }
    }

    if (dim>0) node_matrix = inv(F);

}

template<unsigned int dim, unsigned int spacedim>
FEInternalData *FE_RT0<dim,spacedim>::initialize(const Quadrature<dim> &q, UpdateFlags flags)
{
    FEInternalData *data = new FEInternalData;

    if (flags & update_values)
    {
        arma::mat::fixed<n_raw_functions,dim> raw_values;
        arma::mat::fixed<dim+1,dim> shape_values;
        vector<arma::vec> values;

        data->basis_vectors.resize(q.size());
        values.resize(dim+1);
        for (unsigned int i=0; i<q.size(); i++)
        {
            for (unsigned int j=0; j<n_raw_functions; j++)
                for (unsigned int c=0; c<dim; c++)
                  raw_values(j,c) = basis_value(j, q.point(i), c);

            shape_values = node_matrix * raw_values;

            for (unsigned int j=0; j<dim+1; j++)
                values[j] = trans(shape_values.row(j));

            data->basis_vectors[i] = values;
        }
    }

    if (flags & update_gradients)
    {
        arma::mat::fixed<dim,dim> grad;
        arma::mat::fixed<dim,dim> shape_grads;
        vector<arma::mat> grads;

        data->basis_grad_vectors.resize(q.size());
        grads.resize(dim+1);
        for (unsigned int i=0; i<q.size(); i++)
        {
            for (unsigned int k=0; k<dim+1; k++)
            {
                grad.zeros();
                for (unsigned int l=0; l<n_raw_functions; l++)
                  for (unsigned int c=0; c<dim; c++)
                    grad.col(c) += basis_grad(l, q.point(i), c) * node_matrix(k,l);
                grads[k] = grad;
            }

            data->basis_grad_vectors[i] = grads;
        }
    }

    return data;
}

template<unsigned int dim, unsigned int spacedim> inline
UpdateFlags FE_RT0<dim,spacedim>::update_each(UpdateFlags flags)
{
    UpdateFlags f = flags;

    if (flags & update_values)
        f |= update_jacobians | update_volume_elements;

    if (flags & update_gradients)
        f |= update_jacobians | update_inverse_jacobians | update_volume_elements;

    return f;
}

template<unsigned int dim, unsigned int spacedim> inline
void FE_RT0<dim,spacedim>::fill_fe_values(
        const Quadrature<dim> &q,
        FEInternalData &data,
        FEValuesData<dim,spacedim> &fv_data)
{
    // shape values
    if (fv_data.update_flags & update_values)
    {
        arma::vec::fixed<spacedim> values;
        for (unsigned int i = 0; i < q.size(); i++)
        {
            for (unsigned int k=0; k<dim+1; k++)
            {
                values = fv_data.jacobians[i]*data.basis_vectors[i][k]/fv_data.determinants[i];
                for (unsigned int c=0; c<spacedim; c++)
                  fv_data.shape_values[i][k*spacedim+c] = values(c);
            }
        }
    }

    // shape gradients
    if (fv_data.update_flags & update_gradients)
    {
        arma::mat::fixed<spacedim,spacedim> grads;
        for (unsigned int i = 0; i < q.size(); i++)
        {
            for (unsigned int k=0; k<dim+1; k++)
            {
              grads = fv_data.jacobians[i]*data.basis_grad_vectors[i][k]*fv_data.inverse_jacobians[i]/fv_data.determinants[i];
              for (unsigned int c=0; c<spacedim; c++)
                fv_data.shape_gradients[i][k*spacedim+c] = grads.col(c);
            }
        }
    }
}


template class FE_RT0<1,3>;
template class FE_RT0<2,3>;
template class FE_RT0<3,3>;

