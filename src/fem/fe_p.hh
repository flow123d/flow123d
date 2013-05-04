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
 * @brief Definitions of basic Lagrangean finite elements with polynomial shape functions.
 * @author Jan Stebel
 */

#ifndef FE_P_HH_
#define FE_P_HH_

#include <vector>

#include "system/global_defs.h"
#include "system/system.hh"
#include "fem/finite_element.hh"


/**
 * @brief Space of polynomial functions.
 *
 * This class serves for evaluation of the value and gradient
 * of a polynomial of order @p degree in @p dim variables.
 */
template<unsigned int degree, unsigned int dim>
class PolynomialSpace
{
public:

	/**
	 * @brief Constructor.
	 *
	 * Creates the coefficients of the basis.
	 */
    PolynomialSpace();

    /**
     * @brief Value of the @p i th basis function at point @p p.
     * @param i Number of the basis function.
     * @param p Point at which the function is evaluated.
     */
    const double basis_value(unsigned int i, const arma::vec::fixed<dim> &p) const;

    /**
     * @brief Gradient of the @p i th basis function at point @p p.
     * @param i Number of the basis function.
     * @param p Point at which the function is evaluated.
     */
    const arma::vec::fixed<dim> basis_grad(unsigned int i, const arma::vec::fixed<dim> &p) const;

private:

    /**
     * @brief Coefficients of basis functions.
     *
     * Powers of x, y, z, ... in the i-th basis function are stored
     * in powers[i].
     */
    std::vector<arma::uvec::fixed<dim> > powers;

};



/**
 * @brief Distribution of dofs for polynomial finite elements.
 *
 * The class holds the information on the total number of dofs
 * as well as the number of dofs associated to geometrical entities
 * such as points, lines, triangles and tetrahedra.
 * Moreover, some dofs are grouped to pairs, triples or sextuples
 * which are invariant to rotation/reflection of the element.
 *
 * The coordinates of unit support points are provided.
 * The values at support points uniquely determine the finite
 * element function.
 *
 */
template<unsigned int degree, unsigned int dim>
class DofDistribution
{
public:

	/**
	 * @brief Constructor.
	 *
	 * Initializes all variables.
	 */
    DofDistribution();

    /// Total number of degrees of freedom at one finite element.
    unsigned int number_of_dofs;

    /**
     * @brief Number of single dofs at one geometrical entity of the given
     * dimension (point, line, triangle, tetrahedron).
     */
    unsigned int number_of_single_dofs[dim + 1];

    /**
     * @brief Number of pairs of dofs at one geometrical entity of the given
     * dimension (applicable to lines and triangles).
     */
    unsigned int number_of_pairs[dim + 1];

    /**
     * @brief Number of triples of dofs associated to one triangle.
     */
    unsigned int number_of_triples[dim + 1];

    /**
     * @brief Number of sextuples of dofs associated to one triangle.
     */
    unsigned int number_of_sextuples[dim + 1];

    /**
     * @brief Support points.
     *
     * Support points are points in the reference element where
     * function values determine the dofs. In case of Lagrangean
     * finite elements the dof values are precisely the function
     * values at @p unit_support_points.
     */
    std::vector<arma::vec::fixed<dim> > unit_support_points;

};


/**
 * @brief Conforming Lagrangean finite element on @p dim dimensional simplex.
 *
 * The finite element functions are continuous across the interfaces.
 */
template <unsigned int degree, unsigned int dim, unsigned int spacedim>
class FE_P : public FiniteElement<dim,spacedim>
{
    using FiniteElement<dim,spacedim>::number_of_dofs;
    using FiniteElement<dim,spacedim>::number_of_single_dofs;
    using FiniteElement<dim,spacedim>::number_of_pairs;
    using FiniteElement<dim,spacedim>::number_of_triples;
    using FiniteElement<dim,spacedim>::number_of_sextuples;
    using FiniteElement<dim,spacedim>::unit_support_points;
    using FiniteElement<dim,spacedim>::order;

public:
    /// Constructor.
    FE_P();

    /**
     * @brief Returns the @p ith basis function evaluated at the point @p p.
     * @param i Number of the basis function.
     * @param p Point of evaluation.
     */
    double basis_value(const unsigned int i, const arma::vec::fixed<dim> &p) const;

    /**
     * @brief Returns the gradient of the @p ith basis function at the point @p p.
     * @param i Number of the basis function.
     * @param p Point of evaluation.
     */
    arma::vec::fixed<dim> basis_grad(const unsigned int i, const arma::vec::fixed<dim> &p) const;

    /**
     * @brief The vector variant of basis_value must be implemented but may not be used.
     */
    arma::vec::fixed<dim> basis_vector(const unsigned int i, const arma::vec::fixed<dim> &p) const;

    /**
     * @brief The vector variant of basis_grad must be implemented but may not be used.
     */
    arma::mat::fixed<dim,dim> basis_grad_vector(const unsigned int i, const arma::vec::fixed<dim> &p) const;

private:

    /// The auxiliary polynomial space.
    PolynomialSpace<degree,dim> poly_space;

    /// The auxiliary dof distribution.
    DofDistribution<degree,dim> dof_distribution;
};


/**
 * @brief Discontinuous Lagrangean finite element on @p dim dimensional simplex.
 *
 * No continuity of the finite element functions across the interfaces is
 * imposed.
 */
template <unsigned int degree, unsigned int dim, unsigned int spacedim>
class FE_P_disc : public FiniteElement<dim,spacedim>
{
    using FiniteElement<dim,spacedim>::number_of_dofs;
    using FiniteElement<dim,spacedim>::number_of_single_dofs;
    using FiniteElement<dim,spacedim>::number_of_pairs;
    using FiniteElement<dim,spacedim>::number_of_triples;
    using FiniteElement<dim,spacedim>::number_of_sextuples;
    using FiniteElement<dim,spacedim>::unit_support_points;
    using FiniteElement<dim,spacedim>::order;

public:

    /// Constructor.
    FE_P_disc();

    /**
     * @brief Returns the @p ith basis function evaluated at the point @p p.
     * @param i Number of the basis function.
     * @param p Point of evaluation.
     */
    double basis_value(const unsigned int i, const arma::vec::fixed<dim> &p) const;

    /**
     * @brief Returns the gradient of the @p ith basis function at the point @p p.
     * @param i Number of the basis function.
     * @param p Point of evaluation.
     */
    arma::vec::fixed<dim> basis_grad(const unsigned int i, const arma::vec::fixed<dim> &p) const;

    /**
     * @brief The vector variant of basis_value must be implemented but may not be used.
     */
    arma::vec::fixed<dim> basis_vector(const unsigned int i, const arma::vec::fixed<dim> &p) const;

    /**
     * @brief The vector variant of basis_grad must be implemented but may not be used.
     */
    arma::mat::fixed<dim,dim> basis_grad_vector(const unsigned int i, const arma::vec::fixed<dim> &p) const;

private:

    /// The auxiliary polynomial space.
    PolynomialSpace<degree,dim> poly_space;

    /// The auxiliary dof distribution.
    DofDistribution<degree,dim> dof_distribution;
};



template<unsigned int degree, unsigned int dim>
PolynomialSpace<degree,dim>::PolynomialSpace()
{
// computes powers of all monomials up to given @p degree
// the order is: 1,x,x^2, y, yx,y^2
//
// TODO: - check and possibly rewrite to be more clear (use sum_degree temporary
//       - change order of monomials: 1, x, y, xy, x^2 , y^2 (increasing order)
//       - allow Q polynomials: 1,x, y, xy
//       - can use tensor products

	arma::uvec::fixed<dim> pows;
	pows.zeros();

    unsigned int degree_sum=0;
    unsigned int i_dim;


    while (true) {
        powers.push_back(pows);

        // increment pows
        for(i_dim=0; i_dim < dim; i_dim++) {
            if (degree_sum < degree) {
                pows[i_dim]++;
                degree_sum++;
                break;
            } else {                    // if degree_sum == degree, we find first non empty power, free it, and raise the next one
                degree_sum-=pows[i_dim];
                pows[i_dim]=0;
            }
        }
        if (i_dim == dim) break; // just after pow == (0, 0, .., degree)
    }
}

template<unsigned int degree, unsigned int dim>
const double PolynomialSpace<degree,dim>::basis_value(unsigned int i, const arma::vec::fixed<dim> &p) const
{
    ASSERT(i<=powers.size(), "Index of basis function is out of range.");

    double v = 1;

    for (unsigned int j=0; j<dim; j++)
        v *= pow(p[j], (int) powers[i][j]);

    return v;
}


template<unsigned int degree, unsigned int dim>
const arma::vec::fixed<dim> PolynomialSpace<degree,dim>::basis_grad(unsigned int i, const arma::vec::fixed<dim> &p) const
{
    ASSERT(i<=powers.size(), "Index of basis function is out of range.");

    arma::vec::fixed<dim> grad;

    for (unsigned int j=0; j<dim; j++)
    {
        grad[j] = powers[i][j];
        if (powers[i][j] == 0) continue;

        for (unsigned int k=0; k<dim; k++)
        {
            grad[j] *= pow(p[k], (int) (k==j?powers[i][k]-1:powers[i][k]));
        }
    }
    return grad;
}











template<unsigned int degree, unsigned int dim, unsigned int spacedim>
FE_P<degree,dim,spacedim>::FE_P()
{
    this->init();

    for (int i=0; i<=dim; i++)
    {
        number_of_dofs += dof_distribution.number_of_single_dofs[i]
                         +dof_distribution.number_of_pairs[i]
                         +dof_distribution.number_of_triples[i]
                         +dof_distribution.number_of_sextuples[i];

        number_of_single_dofs[i] = dof_distribution.number_of_single_dofs[i];
        number_of_pairs[i] = dof_distribution.number_of_pairs[i];
        number_of_triples[i] = dof_distribution.number_of_triples[i];
        number_of_sextuples[i] = dof_distribution.number_of_sextuples[i];
    }

    for (int i=0; i<dof_distribution.unit_support_points.size(); i++)
        unit_support_points.push_back(dof_distribution.unit_support_points[i]);

    order = degree;

    this->compute_node_matrix();
}

template<unsigned int degree, unsigned int dim, unsigned int spacedim>
double FE_P<degree,dim,spacedim>::basis_value(const unsigned int i, const arma::vec::fixed<dim> &p) const
{
    ASSERT(i <= number_of_dofs, "Index of basis function is out of range.");
    return poly_space.basis_value(i, p);
}

template<unsigned int degree, unsigned int dim, unsigned int spacedim>
arma::vec::fixed<dim> FE_P<degree,dim,spacedim>::basis_grad(const unsigned int i, const arma::vec::fixed<dim> &p) const
{
    ASSERT(i <= number_of_dofs, "Index of basis function is out of range.");
    return poly_space.basis_grad(i, p);
}

template<unsigned int degree, unsigned int dim, unsigned int spacedim>
arma::vec::fixed<dim> FE_P<degree,dim,spacedim>::basis_vector(const unsigned int i, const arma::vec::fixed<dim> &p) const
{
    ASSERT(false, "basis_vector() may not be called for scalar finite element.");
}

template<unsigned int degree, unsigned int dim, unsigned int spacedim>
arma::mat::fixed<dim,dim> FE_P<degree,dim,spacedim>::basis_grad_vector(const unsigned int i, const arma::vec::fixed<dim> &p) const
{
    ASSERT(false, "basis_grad_vector() may not be called for scalar finite element.");
}















template<unsigned int degree, unsigned int dim, unsigned int spacedim>
FE_P_disc<degree,dim,spacedim>::FE_P_disc()
{
    this->init();

    number_of_dofs += dof_distribution.number_of_dofs;

    number_of_single_dofs[dim] = number_of_dofs;

    for (unsigned int i=0; i<dof_distribution.unit_support_points.size(); i++)
        unit_support_points.push_back(dof_distribution.unit_support_points[i]);

    order = degree;

    this->compute_node_matrix();
}

template<unsigned int degree, unsigned int dim, unsigned int spacedim>
double FE_P_disc<degree,dim,spacedim>::basis_value(const unsigned int i, const arma::vec::fixed<dim> &p) const
{
    ASSERT(i <= number_of_dofs, "Index of basis function is out of range.");
    return poly_space.basis_value(i, p);
}

template<unsigned int degree, unsigned int dim, unsigned int spacedim>
arma::vec::fixed<dim> FE_P_disc<degree,dim,spacedim>::basis_grad(const unsigned int i, const arma::vec::fixed<dim> &p) const
{
    ASSERT(i <= number_of_dofs, "Index of basis function is out of range.");
    return poly_space.basis_grad(i, p);
}

template<unsigned int degree, unsigned int dim, unsigned int spacedim>
arma::vec::fixed<dim> FE_P_disc<degree,dim,spacedim>::basis_vector(const unsigned int i, const arma::vec::fixed<dim> &p) const
{
    ASSERT(false, "basis_vector() may not be called for scalar finite element.");
}

template<unsigned int degree, unsigned int dim, unsigned int spacedim>
arma::mat::fixed<dim,dim> FE_P_disc<degree,dim,spacedim>::basis_grad_vector(const unsigned int i, const arma::vec::fixed<dim> &p) const
{
    ASSERT(false, "basis_grad_vector() may not be called for scalar finite element.");
}












#ifndef DOXYGEN_SHOULD_SKIP_THIS

/****** Template specializations ******/

/*** 1D finite elements ***/

// P0 constant element
template<>
DofDistribution<0,1>::DofDistribution();

// P1 linear element
template<>
DofDistribution<1,1>::DofDistribution();

// P2 quadratic element
template<>
DofDistribution<2,1>::DofDistribution();

/*** 2D finite elements ***/

// P0 constant element
template<>
DofDistribution<0,2>::DofDistribution();


// P1 linear element
template<>
DofDistribution<1,2>::DofDistribution();

// P2 quadratic element
template<>
DofDistribution<2,2>::DofDistribution();



/*** 3D finite elements ***/

// P0 constant element
template<>
DofDistribution<0,3>::DofDistribution();


// P1 linear element
template<>
DofDistribution<1,3>::DofDistribution();


// P2 quadratic element
template<>
DofDistribution<2,3>::DofDistribution();

#endif /* DOXYGEN_SHOULD_SKIP_THIS */




#endif /* FE_P_HH_ */
