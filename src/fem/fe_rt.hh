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
 * @brief Definitions of Raviart-Thomas finite elements.
 * @author Jan Stebel
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

};










template<unsigned int dim, unsigned int spacedim>
FE_RT0<dim,spacedim>::FE_RT0()
{
	arma::vec::fixed<dim> sp;

    this->init();

    number_of_dofs = dim+1;
    number_of_single_dofs[dim] = dim+1;

    sp.fill(1./max(1.,(double)dim));
    generalized_support_points.push_back(sp);
    for (int i=0; i<dim; i++)
    {
        sp.fill(1./max(1.,(double)dim));
        sp[i] = 0;
        generalized_support_points.push_back(sp);
    }

    order = 1;

    is_scalar_fe = false;

    compute_node_matrix();
}

template<unsigned int dim, unsigned int spacedim>
double FE_RT0<dim,spacedim>::basis_value(const unsigned int i, const arma::vec::fixed<dim> &p) const
{
    ASSERT(false, "basis_value() may not be called for vectorial finite element.");
}

template<unsigned int dim, unsigned int spacedim>
arma::vec::fixed<dim> FE_RT0<dim,spacedim>::basis_grad(const unsigned int i, const arma::vec::fixed<dim> &p) const
{
    ASSERT(false, "basis_grad() may not be called for vectorial finite element.");
}

template<unsigned int dim, unsigned int spacedim>
arma::vec::fixed<dim> FE_RT0<dim,spacedim>::basis_vector(const unsigned int i, const arma::vec::fixed<dim> &p) const
{
    ASSERT(i<(dim+1)*dim, "Index of basis function is out of range.");

    arma::vec::fixed<dim> v;
    unsigned int comp = i/(dim+1);
    unsigned int offset = i%(dim+1);
    
    v.zeros();

    if (offset < dim)
    {
        v[comp] = p[offset];
    }
    else
    {
        v[comp] = 1;
    }

    return v;
}

template<unsigned int dim, unsigned int spacedim>
arma::mat::fixed<dim,dim> FE_RT0<dim,spacedim>::basis_grad_vector(const unsigned int i, const arma::vec::fixed<dim> &p) const
{
    ASSERT(i<(dim+1)*dim, "Index of basis function is out of range.");

    arma::mat::fixed<dim,dim> m;
    unsigned int comp = i/(dim+1);
    unsigned int offset = i%(dim+1);

    m.zeros();

    if (offset < dim)
    {
        m(comp,offset) = 1;
    }

    return m;
}

template<unsigned int dim, unsigned int spacedim>
void FE_RT0<dim,spacedim>::compute_node_matrix()
{
	arma::mat::fixed<dim*(dim+1),dim+1> F;
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

    for (int i=0; i<dim*(dim+1); i++)
    {
        /*
         * For the 0-th side we have:
         *
         * $$ |\Gamma_0| = \frac{\sqrt{dim}}{dim-1}, $$
         * $$ n_0 = \frac{\sqrt{dim}}{dim}(1,\ldots,1). $$
         *
         * If dim==1, we set |\Gamma_0| = 1.
         *
         */

        r = basis_vector(i,generalized_support_points[0]);
        if (dim == 1)
        {
            F(i,0) = sum(r);
        }
        else
        {
            F(i,0) = sum(r)/(1.*dim-1);
        }

        /*
         * For the remaining sides:
         *
         * $$ |\Gamma_j| = \sqrt(dim-1), $$
         * $$ (n_j)_l = -\delta_{jl}. $$
         *
         * If dim==1, we set |\Gamma_1| = 1.
         *
         */

        if (dim == 1)
        {
            r = basis_vector(i,generalized_support_points[1]);
            F(i,1) = -r(0);
        }
        else
        {
            for (int j=1; j<dim+1; j++)
            {
                r = basis_vector(i,generalized_support_points[j]);
                F(i,j) = -r(j-1)/(1.*dim-1);
            }
        }
    }

    if (dim>0) node_matrix = pinv(F);
}

template<unsigned int dim, unsigned int spacedim>
FEInternalData *FE_RT0<dim,spacedim>::initialize(const Quadrature<dim> &q, UpdateFlags flags)
{
    FEInternalData *data = new FEInternalData;

    if (flags & update_values)
    {
    	arma::mat::fixed<dim*(dim+1),dim> raw_values;
    	arma::mat::fixed<dim+1,dim> shape_values;
        vector<arma::vec> values;

        data->basis_vectors.resize(q.size());
        values.resize(dim+1);
        for (int i=0; i<q.size(); i++)
        {
            for (int j=0; j<dim*(dim+1); j++)
                raw_values.row(j) = trans(basis_vector(j, q.point(i)));

            shape_values = node_matrix * raw_values;

            for (int j=0; j<dim+1; j++)
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
        for (int i=0; i<q.size(); i++)
        {
            for (int k=0; k<dim+1; k++)
            {
                grad.zeros();
                for (int l=0; l<dim*(dim+1); l++)
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
        for (int i = 0; i < q.size(); i++)
        {
            for (int k=0; k<dim+1; k++)
                vectors[k] = fv_data.jacobians[i]*data.basis_vectors[i][k]/fv_data.determinants[i];

            fv_data.shape_vectors[i] = vectors;
        }
    }

    // shape gradients
    if (fv_data.update_flags & update_gradients)
    {
        vector<arma::mat::fixed<spacedim,spacedim> > grads;
        grads.resize(dim+1);
        for (int i = 0; i < q.size(); i++)
        {
            for (int k=0; k<dim+1; k++)
                grads[k] = fv_data.jacobians[i]*data.basis_grad_vectors[i][k]*fv_data.inverse_jacobians[i]/fv_data.determinants[i];

            fv_data.shape_grad_vectors[i] = grads;
        }
    }
}













#endif /* FE_RT_HH_ */
