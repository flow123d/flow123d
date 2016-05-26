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
 * @file    fe_rt.hh
 * @brief   Definitions of Raviart-Thomas finite elements.
 * @author  Jan Stebel
 */

#ifndef FE_RT_HH_
#define FE_RT_HH_

#include "fem/finite_element.hh"



/**
 * @brief Raviart-Thomas element of order 0.
 *
 * The lowest order Raviart-Thomas finite element with linear basis functions
 * and continuous normal components across element sides.
 */
template <unsigned int dim, unsigned int spacedim>
class FE_RT0 : public FiniteElement<dim,spacedim>
{
    using FiniteElement<dim,spacedim>::number_of_dofs;
    using FiniteElement<dim,spacedim>::number_of_single_dofs;
    using FiniteElement<dim,spacedim>::number_of_pairs;
    using FiniteElement<dim,spacedim>::number_of_triples;
    using FiniteElement<dim,spacedim>::number_of_sextuples;
    using FiniteElement<dim,spacedim>::generalized_support_points;
    using FiniteElement<dim,spacedim>::order;
    using FiniteElement<dim,spacedim>::is_scalar_fe;
    using FiniteElement<dim,spacedim>::node_matrix;

public:

    /// Number of raw basis functions.
    static const unsigned int n_raw_functions = dim+1;

    /**
     * @brief Constructor.
     */
    FE_RT0();

    /**
     * @brief The scalar variant of basis_vector must be implemented but may not be used.
     */
    double basis_value(const unsigned int i, const arma::vec::fixed<dim> &p) const;

    /**
     * @brief The scalar variant of basis_grad_vector must be implemented but may not be used.
     */
    arma::vec::fixed<dim> basis_grad(const unsigned int i, const arma::vec::fixed<dim> &p) const;

    /**
     * @brief Returns the @p ith basis function evaluated at the point @p p.
     * @param i Number of the basis function.
     * @param p Point of evaluation.
     */
    arma::vec::fixed<dim> basis_vector(const unsigned int i, const arma::vec::fixed<dim> &p) const;

    /**
     * @brief Returns the gradient of the @p ith basis function at the point @p p.
     * @param i Number of the basis function.
     * @param p Point of evaluation.
     */
    arma::mat::fixed<dim,dim> basis_grad_vector(const unsigned int i, const arma::vec::fixed<dim> &p) const;

    /**
     * @brief Computes the conversion matrix from internal basis to shape functions.
     *
     * Initializes the @p node_matrix for computing the coefficients
     * of the raw basis functions from values at support points.
     * Since Raviart-Thomas element is not Lagrangean, the method has
     * to be reimplemented.
     */
    void compute_node_matrix();

    /**
     * @brief Calculates the data on the reference cell.
     *
     * @param q Quadrature.
     * @param flags Flags that indicate what quantities should be calculated.
     */
    FEInternalData *initialize(const Quadrature<dim> &q, UpdateFlags flags);

    /**
     * @brief Decides which additional quantities have to be computed
     * for each cell.
     */
    UpdateFlags update_each(UpdateFlags flags);

    /**
     * @brief Computes the shape function values and gradients on the actual cell
     * and fills the FEValues structure.
     *
     * @param q Quadrature.
     * @param data The precomputed finite element data on the reference cell.
     * @param fv_data The data to be computed.
     */
    virtual void fill_fe_values(const Quadrature<dim> &q,
            FEInternalData &data,
            FEValuesData<dim,spacedim> &fv_data);

    /**
     * @brief Destructor.
     */
    virtual ~FE_RT0() {};

};










template<unsigned int dim, unsigned int spacedim>
FE_RT0<dim,spacedim>::FE_RT0()
{
	arma::vec::fixed<dim> sp;

    this->init();

    number_of_dofs = dim+1;
    number_of_single_dofs[dim] = dim+1;

    for (unsigned int sid=0; sid<RefElement<dim>::n_sides; ++sid)
    {
    	sp.fill(0);
    	for (unsigned int i=0; i<RefElement<dim>::n_nodes_per_side; ++i)
    		sp += RefElement<dim>::node_coords(RefElement<dim>::side_nodes[sid][i]);
    	sp /= RefElement<dim>::n_nodes_per_side;
    	generalized_support_points.push_back(sp);
    }

    order = 1;

    is_scalar_fe = false;

    compute_node_matrix();
}

template<unsigned int dim, unsigned int spacedim>
double FE_RT0<dim,spacedim>::basis_value(const unsigned int i, const arma::vec::fixed<dim> &p) const
{
	OLD_ASSERT(false, "basis_value() may not be called for vectorial finite element.");

    return 0.0;
}

template<unsigned int dim, unsigned int spacedim>
arma::vec::fixed<dim> FE_RT0<dim,spacedim>::basis_grad(const unsigned int i, const arma::vec::fixed<dim> &p) const
{
	OLD_ASSERT(false, "basis_grad() may not be called for vectorial finite element.");
    return arma::vec::fixed<dim>();
}

template<unsigned int dim, unsigned int spacedim>
arma::vec::fixed<dim> FE_RT0<dim,spacedim>::basis_vector(const unsigned int i, const arma::vec::fixed<dim> &p) const
{
	OLD_ASSERT(i<n_raw_functions, "Index of basis function is out of range.");

    arma::vec::fixed<dim> v(p);
    
    if (i > 0)
    	v[i-1] -= 1;

    return v;
}

template<unsigned int dim, unsigned int spacedim>
arma::mat::fixed<dim,dim> FE_RT0<dim,spacedim>::basis_grad_vector(const unsigned int i, const arma::vec::fixed<dim> &p) const
{
    OLD_ASSERT(i<n_raw_functions, "Index of basis function is out of range.");

    return arma::eye(dim,dim);
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
        	r = basis_vector(i,generalized_support_points[j]);
        	F(i,j) = dot(r,RefElement<dim>::normal_vector(j))*RefElement<dim>::side_measure(j);
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
                raw_values.row(j) = trans(basis_vector(j, q.point(i)));

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
                    grad += basis_grad_vector(l, q.point(i)) * node_matrix(k,l);
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
        vector<arma::vec::fixed<spacedim> > vectors;
        vectors.resize(dim+1);
        for (unsigned int i = 0; i < q.size(); i++)
        {
            for (unsigned int k=0; k<dim+1; k++)
                vectors[k] = fv_data.jacobians[i]*data.basis_vectors[i][k]/fv_data.determinants[i];

            fv_data.shape_vectors[i] = vectors;
        }
    }

    // shape gradients
    if (fv_data.update_flags & update_gradients)
    {
        vector<arma::mat::fixed<spacedim,spacedim> > grads;
        grads.resize(dim+1);
        for (unsigned int i = 0; i < q.size(); i++)
        {
            for (unsigned int k=0; k<dim+1; k++)
                grads[k] = fv_data.jacobians[i]*data.basis_grad_vectors[i][k]*fv_data.inverse_jacobians[i]/fv_data.determinants[i];

            fv_data.shape_grad_vectors[i] = grads;
        }
    }
}













#endif /* FE_RT_HH_ */
