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
 * $Id: quadrature.hh 1352 2011-09-23 14:14:47Z jan.stebel $
 * $Revision: 1352 $
 * $LastChangedBy: jan.stebel $
 * $LastChangedDate: 2011-09-23 16:14:47 +0200 (Fri, 23 Sep 2011) $
 *
 * @file
 * @brief Definitions of basic Lagrangean finite elements with polynomial shape functions.
 * @author Jan Stebel
 */

#ifndef FE_P_HH_
#define FE_P_HH_

#include "fem/finite_element.hh"

using namespace arma;


/**
 * Conforming Lagrangean finite element on @p dim dimensional simplex.
 * The finite element functions are continuous across the interfaces.
 */
template <unsigned int dim, unsigned int degree>
class FE_P : public FiniteElement<dim>
{
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
};

/**
 * Discontinuous Lagrangean finite element on @p dim dimensional simplex.
 * No continuity of the finite element functions across the interfaces is
 * imposed.
 */
template <unsigned int dim, unsigned int degree>
class FE_P_disc : public FiniteElement<dim>
{
private:
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
};



/****** Template specializations ******/

/*** 1D finite elements ***/

// P0 constant element
template<>
FE_P<1,0>::FE_P()
{
    number_of_dofs = 1;

    number_of_single_dofs[1] = 1;

    unit_support_points.push_back(zeros<vec>(1));

    compute_node_matrix();
}

template<>
double FE_P<1,0>::basis_value(const unsigned int i, const vec::fixed<1> &p) const
{
    ASSERT(i==0, "Index of shape function is out of range.");
    return 1;
}

template<>
vec::fixed<1> FE_P<1,0>::basis_grad(const unsigned int i, const vec::fixed<1> &p) const
{
    ASSERT(i==0, "Index of shape function is out of range.");
    return zeros<vec>(1);
}


// P1 linear element
template<>
FE_P<1,1>::FE_P()
{
    number_of_dofs = 2;

    number_of_single_dofs[0] = 2;

    unit_support_points.push_back(vec::fixed<1>("0"));
    unit_support_points.push_back(vec::fixed<1>("1"));
    compute_node_matrix();
}

template<>
double FE_P<1,1>::basis_value(const unsigned int i, const vec::fixed<1> &p) const
{
    ASSERT(i <= 1, "Index of shape function is out of range.");
    switch (i)
    {
    case 0:
        return 1;
        break;
    case 1:
        return p(0);
        break;
    }
}

template<>
vec::fixed<1> FE_P<1,1>::basis_grad(const unsigned int i, const vec::fixed<1> &p) const
{
    ASSERT(i <= 1, "Index of shape function is out of range.");
    static const vec grad = ones<vec>(1);
    switch (i)
    {
    case 0:
        return vec::fixed<1>("0");
        break;
    case 1:
        return grad;
        break;
    }
}


// P1 linear discontinuous element
template<>
FE_P_disc<1,1>::FE_P_disc()
{
    number_of_dofs = 2;

    number_of_single_dofs[1] = 2;

    unit_support_points.push_back(vec::fixed<1>("0"));
    unit_support_points.push_back(vec::fixed<1>("1"));
    compute_node_matrix();
}

template<>
double FE_P_disc<1,1>::basis_value(const unsigned int i, const vec::fixed<1> &p) const
{
    ASSERT(i <= 1, "Index of shape function is out of range.");
    switch (i)
    {
    case 0:
        return 1;
        break;
    case 1:
        return p(0);
        break;
    }
}

template<>
vec::fixed<1> FE_P_disc<1,1>::basis_grad(const unsigned int i, const vec::fixed<1> &p) const
{
    ASSERT(i <= 1, "Index of shape function is out of range.");
    static const vec grad = ones<vec>(1);
    switch (i)
    {
    case 0:
        return vec::fixed<1>("0");
        break;
    case 1:
        return grad;
        break;
    }
}


/*** 2D finite elements ***/

// P0 constant element
template<>
FE_P<2,0>::FE_P()
{
    number_of_dofs = 1;

    number_of_single_dofs[2] = 1;

    unit_support_points.push_back(vec2("0 0"));
    compute_node_matrix();
}

template<>
double FE_P<2,0>::basis_value(const unsigned int i, const vec2 &p) const
{
    ASSERT(i==0, "Index of shape function is out of range.");
    return 1;
}

template<>
vec2 FE_P<2,0>::basis_grad(const unsigned int i, const vec2 &p) const
{
    ASSERT(i==0, "Index of shape function is out of range.");
    return zeros<vec>(2);
}


// P1 linear element
template<>
FE_P<2,1>::FE_P()
{
    number_of_dofs = 3;

    number_of_single_dofs[0] = 3;

    unit_support_points.push_back(vec2("0 0"));
    unit_support_points.push_back(vec2("1 0"));
    unit_support_points.push_back(vec2("0 1"));
    compute_node_matrix();
}

template<>
double FE_P<2,1>::basis_value(const unsigned int i, const vec2 &p) const
{
    ASSERT(i <= 2, "Index of shape function is out of range.");
    switch (i)
    {
    case 0:
        return 1;
        break;
    case 1:
        return p(0);
        break;
    case 2:
        return p(1);
        break;
    }
}

template<>
vec2 FE_P<2,1>::basis_grad(const unsigned int i, const vec2 &p) const
{
    ASSERT(i <= 2, "Index of shape function is out of range.");
    switch (i)
    {
    case 0:
        return "0 0";
        break;
    case 1:
        return "1 0";
        break;
    case 2:
        return "0 1";
        break;
    }
}


// P1 linear discontinuous element
template<>
FE_P_disc<2,1>::FE_P_disc()
{
    number_of_dofs = 3;

    number_of_single_dofs[2] = 3;

    unit_support_points.push_back(vec2("0 0"));
    unit_support_points.push_back(vec2("1 0"));
    unit_support_points.push_back(vec2("0 1"));
    compute_node_matrix();
}

template<>
double FE_P_disc<2,1>::basis_value(const unsigned int i, const vec2 &p) const
{
    ASSERT(i <= 2, "Index of shape function is out of range.");
    switch (i)
    {
    case 0:
        return 1;
        break;
    case 1:
        return p(0);
        break;
    case 2:
        return p(1);
        break;
    }
}

template<>
vec2 FE_P_disc<2,1>::basis_grad(const unsigned int i, const vec2 &p) const
{
    ASSERT(i <= 2, "Index of shape function is out of range.");
    switch (i)
    {
    case 0:
        return "0 0";
        break;
    case 1:
        return "1 0";
        break;
    case 2:
        return "0 1";
        break;
    }
}


/*** 3D finite elements ***/

// P0 constant element
template<>
FE_P<3,0>::FE_P()
{
    number_of_dofs = 1;

    number_of_single_dofs[3] = 1;

    unit_support_points.push_back(vec3("0 0 0"));
    compute_node_matrix();
}

template<>
double FE_P<3,0>::basis_value(const unsigned int i, const vec3 &p) const
{
    ASSERT(i==0, "Index of shape function is out of range.");
    return 1;
}

template<>
vec3 FE_P<3,0>::basis_grad(const unsigned int i, const vec3 &p) const
{
    ASSERT(i==0, "Index of shape function is out of range.");
    return zeros<vec>(3);
}


// P1 linear element
template<>
FE_P<3,1>::FE_P()
{
    number_of_dofs = 4;

    number_of_single_dofs[0] = 4;

    unit_support_points.push_back(vec3("0 0 0"));
    unit_support_points.push_back(vec3("1 0 0"));
    unit_support_points.push_back(vec3("0 1 0"));
    unit_support_points.push_back(vec3("0 0 1"));
    compute_node_matrix();
}

template<>
double FE_P<3,1>::basis_value(const unsigned int i, const vec3 &p) const
{
    ASSERT(i <= 3, "Index of shape function is out of range.");
    switch (i)
    {
    case 0:
        return 1;
        break;
    case 1:
        return p(0);
        break;
    case 2:
        return p(1);
        break;
    case 3:
        return p(2);
        break;
    }
}

template<>
vec3 FE_P<3,1>::basis_grad(const unsigned int i, const vec3 &p) const
{
    ASSERT(i <= 3, "Index of shape function is out of range.");
    switch (i)
    {
    case 0:
        return "0 0 0";
        break;
    case 1:
        return "1 0 0";
        break;
    case 2:
        return "0 1 0";
        break;
    case 3:
        return "0 0 1";
        break;
    }
}


// P1 linear element
template<>
FE_P_disc<3,1>::FE_P_disc()
{
    number_of_dofs = 4;

    number_of_single_dofs[3] = 4;

    unit_support_points.push_back(vec3("0 0 0"));
    unit_support_points.push_back(vec3("1 0 0"));
    unit_support_points.push_back(vec3("0 1 0"));
    unit_support_points.push_back(vec3("0 0 1"));
    compute_node_matrix();
}

template<>
double FE_P_disc<3,1>::basis_value(const unsigned int i, const vec3 &p) const
{
    ASSERT(i <= 3, "Index of shape function is out of range.");
    switch (i)
    {
    case 0:
        return 1;
        break;
    case 1:
        return p(0);
        break;
    case 2:
        return p(1);
        break;
    case 3:
        return p(2);
        break;
    }
}

template<>
vec3 FE_P_disc<3,1>::basis_grad(const unsigned int i, const vec3 &p) const
{
    ASSERT(i <= 3, "Index of shape function is out of range.");
    switch (i)
    {
    case 0:
        return "0 0 0";
        break;
    case 1:
        return "1 0 0";
        break;
    case 2:
        return "0 1 0";
        break;
    case 3:
        return "0 0 1";
        break;
    }
}



#endif /* FE_P_HH_ */
