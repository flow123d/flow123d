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

#include "fem/finite_element.hh"

using namespace arma;

template<unsigned int degree, unsigned int dim>
class PolynomialSpace
{
public:

    PolynomialSpace();

    const double basis_value(unsigned int i, const vec::fixed<dim> &p) const;

    const vec::fixed<dim> basis_grad(unsigned int i, const vec::fixed<dim> &p) const;

private:

    /**
     * Powers of x, y, z, ... in the i-th basis function are stored
     * in powers[i].
     */
    vector<uvec::fixed<dim> > powers;

};


template<unsigned int degree, unsigned int dim>
class DofDistribution
{
public:

    DofDistribution();

    /**
     * Total number of degrees of freedom at one finite element.
     */
    unsigned int number_of_dofs;

    /**
     * Number of single dofs at one geometrical entity of the given
     * dimension (point, line, triangle, tetrahedron).
     */
    unsigned int number_of_single_dofs[dim + 1];

    /**
     * Number of pairs of dofs at one geometrical entity of the given
     * dimension (applicable to lines and triangles).
     */
    unsigned int number_of_pairs[dim + 1];

    /**
     * Number of triples of dofs associated to one triangle.
     */
    unsigned int number_of_triples[dim + 1];

    /**
     * Number of sextuples of dofs associated to one triangle.
     */
    unsigned int number_of_sextuples[dim + 1];

    /**
     * Support points are points in the reference element where
     * function values determine the dofs. In case of Lagrangean
     * finite elements the dof values are precisely the function
     * values at @p unit_support_points.
     *
     */
    vector<vec::fixed<dim> > unit_support_points;

};


/**
 * Conforming Lagrangean finite element on @p dim dimensional simplex.
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
    /**
     * Constructor.
     */
    FE_P();

    /**
     * Returns the @p ith basis function evaluated at the point @p p.
     */
    double basis_value(const unsigned int i, const vec::fixed<dim> &p) const;

    /**
     * Returns the gradient of the @p ith basis function at the point @p p.
     */
    vec::fixed<dim> basis_grad(const unsigned int i, const vec::fixed<dim> &p) const;

    /**
     * The vector variant of basis_value must be implemented but may not be used.
     */
    vec::fixed<dim> basis_vector(const unsigned int i, const vec::fixed<dim> &p) const;

    /**
     * The vector variant of basis_grad must be implemented but may not be used.
     */
    mat::fixed<dim,dim> basis_grad_vector(const unsigned int i, const vec::fixed<dim> &p) const;

private:

    PolynomialSpace<degree,dim> poly_space;
    DofDistribution<degree,dim> dof_distribution;
};


/**
 * Discontinuous Lagrangean finite element on @p dim dimensional simplex.
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
    /**
     * Constructor.
     */
    FE_P_disc();

    /**
     * Returns the @p ith basis function evaluated at the point @p p.
     */
    double basis_value(const unsigned int i, const vec::fixed<dim> &p) const;

    /**
     * Returns the gradient of the @p ith basis function at the point @p p.
     */
    vec::fixed<dim> basis_grad(const unsigned int i, const vec::fixed<dim> &p) const;

    /**
     * The vector variant of basis_value must be implemented but may not be used.
     */
    vec::fixed<dim> basis_vector(const unsigned int i, const vec::fixed<dim> &p) const;

    /**
     * The vector variant of basis_grad must be implemented but may not be used.
     */
    mat::fixed<dim,dim> basis_grad_vector(const unsigned int i, const vec::fixed<dim> &p) const;

private:

    PolynomialSpace<degree,dim> poly_space;
    DofDistribution<degree,dim> dof_distribution;
};



template<unsigned int degree, unsigned int dim>
PolynomialSpace<degree,dim>::PolynomialSpace()
{
    uvec::fixed<dim> pows;
    int i;

    pows.zeros();
    i = 0;

    while (true)
    {
        powers.push_back(pows);

        if (sum(pows) == degree)
        {
            while (pows[i] == 0) i++;
            pows[i] = 0;
            if (i == dim-1) break;
            pows[i+1]++;
            i = 0;
        }
        else
        {
            pows[i]++;
        }
    }
}

template<unsigned int degree, unsigned int dim>
const double PolynomialSpace<degree,dim>::basis_value(unsigned int i, const vec::fixed<dim> &p) const
{
    ASSERT(i<=powers.size(), "Index of basis function is out of range.");

    double v = 1;

    for (int j=0; j<dim; j++)
        v *= pow(p[j], powers[i][j]);

    return v;
}


template<unsigned int degree, unsigned int dim>
const vec::fixed<dim> PolynomialSpace<degree,dim>::basis_grad(unsigned int i, const vec::fixed<dim> &p) const
{
    ASSERT(i<=powers.size(), "Index of basis function is out of range.");

    vec::fixed<dim> grad;

    for (int j=0; j<dim; j++)
    {
        grad[j] = 1;
        for (int k=0; k<dim; k++)
        {
            if (powers[i][j] == 0)
            {
                grad[j] = 0;
                continue;
            }
            grad[j] *= pow(p[k], (int)powers[i][j]-1);
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
double FE_P<degree,dim,spacedim>::basis_value(const unsigned int i, const vec::fixed<dim> &p) const
{
    ASSERT(i <= number_of_dofs, "Index of basis function is out of range.");
    return poly_space.basis_value(i, p);
}

template<unsigned int degree, unsigned int dim, unsigned int spacedim>
vec::fixed<dim> FE_P<degree,dim,spacedim>::basis_grad(const unsigned int i, const vec::fixed<dim> &p) const
{
    ASSERT(i <= number_of_dofs, "Index of basis function is out of range.");
    return poly_space.basis_grad(i, p);
}

template<unsigned int degree, unsigned int dim, unsigned int spacedim>
vec::fixed<dim> FE_P<degree,dim,spacedim>::basis_vector(const unsigned int i, const vec::fixed<dim> &p) const
{
    ASSERT(false, "basis_vector() may not be called for scalar finite element.");
}

template<unsigned int degree, unsigned int dim, unsigned int spacedim>
mat::fixed<dim,dim> FE_P<degree,dim,spacedim>::basis_grad_vector(const unsigned int i, const vec::fixed<dim> &p) const
{
    ASSERT(false, "basis_grad_vector() may not be called for scalar finite element.");
}















template<unsigned int degree, unsigned int dim, unsigned int spacedim>
FE_P_disc<degree,dim,spacedim>::FE_P_disc()
{
    this->init();

    number_of_dofs += dof_distribution.number_of_dofs;

    number_of_single_dofs[dim] = number_of_dofs;

    for (int i=0; i<dof_distribution.unit_support_points.size(); i++)
        unit_support_points.push_back(dof_distribution.unit_support_points[i]);

    order = degree;

    this->compute_node_matrix();
}

template<unsigned int degree, unsigned int dim, unsigned int spacedim>
double FE_P_disc<degree,dim,spacedim>::basis_value(const unsigned int i, const vec::fixed<dim> &p) const
{
    ASSERT(i <= number_of_dofs, "Index of basis function is out of range.");
    return poly_space.basis_value(i, p);
}

template<unsigned int degree, unsigned int dim, unsigned int spacedim>
vec::fixed<dim> FE_P_disc<degree,dim,spacedim>::basis_grad(const unsigned int i, const vec::fixed<dim> &p) const
{
    ASSERT(i <= number_of_dofs, "Index of basis function is out of range.");
    return poly_space.basis_grad(i, p);
}

template<unsigned int degree, unsigned int dim, unsigned int spacedim>
vec::fixed<dim> FE_P_disc<degree,dim,spacedim>::basis_vector(const unsigned int i, const vec::fixed<dim> &p) const
{
    ASSERT(false, "basis_vector() may not be called for scalar finite element.");
}

template<unsigned int degree, unsigned int dim, unsigned int spacedim>
mat::fixed<dim,dim> FE_P_disc<degree,dim,spacedim>::basis_grad_vector(const unsigned int i, const vec::fixed<dim> &p) const
{
    ASSERT(false, "basis_grad_vector() may not be called for scalar finite element.");
}














/****** Template specializations ******/

/*** 1D finite elements ***/

// P0 constant element
template<>
DofDistribution<0,1>::DofDistribution()
{
    number_of_dofs = 1;

    number_of_single_dofs[1] = 1;

    unit_support_points.push_back(zeros<vec>(1));
}

// P1 linear element
template<>
DofDistribution<1,1>::DofDistribution()
{
    number_of_dofs = 2;

    number_of_single_dofs[0] = 2;

    unit_support_points.push_back(vec::fixed<1>("0"));
    unit_support_points.push_back(vec::fixed<1>("1"));
}

/*** 2D finite elements ***/

// P0 constant element
template<>
DofDistribution<0,2>::DofDistribution()
{
    number_of_dofs = 1;

    number_of_single_dofs[2] = 1;

    unit_support_points.push_back(vec2("0 0"));
}


// P1 linear element
template<>
DofDistribution<1,2>::DofDistribution()
{
    number_of_dofs = 3;

    number_of_single_dofs[0] = 3;

    unit_support_points.push_back(vec2("0 0"));
    unit_support_points.push_back(vec2("1 0"));
    unit_support_points.push_back(vec2("0 1"));
}



/*** 3D finite elements ***/

// P0 constant element
template<>
DofDistribution<0,3>::DofDistribution()
{
    number_of_dofs = 1;

    number_of_single_dofs[3] = 1;

    unit_support_points.push_back(vec3("0 0 0"));
}


// P1 linear element
template<>
DofDistribution<1,3>::DofDistribution()
{
    number_of_dofs = 4;

    number_of_single_dofs[0] = 4;

    unit_support_points.push_back(vec3("0 0 0"));
    unit_support_points.push_back(vec3("1 0 0"));
    unit_support_points.push_back(vec3("0 1 0"));
    unit_support_points.push_back(vec3("0 0 1"));
}





#endif /* FE_P_HH_ */
