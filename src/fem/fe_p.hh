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


template <unsigned int dim, unsigned int degree>
class FE_P : FiniteElement<dim>
{
private:
    static const bool _dof_continuity[];
    static const FE_dof_object _dof_objs[];
    static const unsigned int _dof_obj_ids[];
public:
    FE_P() : FiniteElement<dim>((dim==1)?(degree+1):((dim==2)?(degree+1)*(degree+2)/2:(degree+1)*(degree+2)*(degree+3)/6), _dof_continuity, _dof_objs, _dof_obj_ids) {};

    double shape_value(const unsigned int i, const vec::fixed<dim> &p) const;
    vec::fixed<dim> shape_grad(const unsigned int i, const vec::fixed<dim> &p) const;
};



/****** Template specializations ******/

/*** 1D finite elements ***/

// P0 constant element
template<>
const bool FE_P<1,0>::_dof_continuity[] = { false };

template<>
const FE_dof_object FE_P<1,0>::_dof_objs[] = { FE_OBJECT_LINE };

template<>
const unsigned int FE_P<1,0>::_dof_obj_ids[] = { 0 };

template<>
double FE_P<1,0>::shape_value(const unsigned int i, const vec::fixed<1> &p) const
{
    ASSERT(i==0, "Index of shape function is out of range.");
    return 1;
}

template<>
vec::fixed<1> FE_P<1,0>::shape_grad(const unsigned int i, const vec::fixed<1> &p) const
{
    ASSERT(i==0, "Index of shape function is out of range.");
    return zeros<vec>(1);
}


// P1 linear element
template<>
const bool FE_P<1,1>::_dof_continuity[] = { true, true };

template<>
const FE_dof_object FE_P<1,1>::_dof_objs[] = { FE_OBJECT_POINT, FE_OBJECT_POINT };

template<>
const unsigned int FE_P<1,1>::_dof_obj_ids[] = { 0, 1 };

template<>
double FE_P<1,1>::shape_value(const unsigned int i, const vec::fixed<1> &p) const
{
    ASSERT(i <= 1, "Index of shape function is out of range.");
    switch (i)
    {
    case 0:
        return 1-p(0);
        break;
    case 1:
        return p(0);
        break;
    }
}

template<>
vec::fixed<1> FE_P<1,1>::shape_grad(const unsigned int i, const vec::fixed<1> &p) const
{
    ASSERT(i <= 1, "Index of shape function is out of range.");
    static const vec grad = ones<vec>(1);
    switch (i)
    {
    case 0:
        return -grad;
        break;
    case 1:
        return grad;
        break;
    }
}



/*** 2D finite elements ***/

// P0 constant element
template<>
const bool FE_P<2,0>::_dof_continuity[] = { false };

template<>
const FE_dof_object FE_P<2,0>::_dof_objs[] = { FE_OBJECT_TRIANGLE };

template<>
const unsigned int FE_P<2,0>::_dof_obj_ids[] = { 0 };

template<>
double FE_P<2,0>::shape_value(const unsigned int i, const vec2 &p) const
{
    ASSERT(i==0, "Index of shape function is out of range.");
    return 1;
}

template<>
vec2 FE_P<2,0>::shape_grad(const unsigned int i, const vec2 &p) const
{
    ASSERT(i==0, "Index of shape function is out of range.");
    return zeros<vec>(2);
}


// P1 linear element
template<>
const bool FE_P<2,1>::_dof_continuity[] = { true, true, true };

template<>
const FE_dof_object FE_P<2,1>::_dof_objs[] = { FE_OBJECT_POINT, FE_OBJECT_POINT, FE_OBJECT_POINT };

template<>
const unsigned int FE_P<2,1>::_dof_obj_ids[] = { 0, 1, 2 };

template<>
double FE_P<2,1>::shape_value(const unsigned int i, const vec2 &p) const
{
    ASSERT(i <= 2, "Index of shape function is out of range.");
    switch (i)
    {
    case 0:
        return 1-p(0)-p(1);
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
vec2 FE_P<2,1>::shape_grad(const unsigned int i, const vec2 &p) const
{
    ASSERT(i <= 2, "Index of shape function is out of range.");
    switch (i)
    {
    case 0:
        return "-1; -1";
        break;
    case 1:
        return "1; 0";
        break;
    case 2:
        return "0; 1";
        break;
    }
}



/*** 3D finite elements ***/

// P0 constant element
template<>
const bool FE_P<3,0>::_dof_continuity[] = { false };

template<>
const FE_dof_object FE_P<3,0>::_dof_objs[] = { FE_OBJECT_TETRAHEDRON };

template<>
const unsigned int FE_P<3,0>::_dof_obj_ids[] = { 0 };

template<>
double FE_P<3,0>::shape_value(const unsigned int i, const vec3 &p) const
{
    ASSERT(i==0, "Index of shape function is out of range.");
    return 1;
}

template<>
vec3 FE_P<3,0>::shape_grad(const unsigned int i, const vec3 &p) const
{
    ASSERT(i==0, "Index of shape function is out of range.");
    return zeros<vec>(3);
}


// P1 linear element
template<>
const bool FE_P<3,1>::_dof_continuity[] = { true, true, true, true };

template<>
const FE_dof_object FE_P<3,1>::_dof_objs[] = { FE_OBJECT_POINT, FE_OBJECT_POINT, FE_OBJECT_POINT, FE_OBJECT_POINT };

template<>
const unsigned int FE_P<3,1>::_dof_obj_ids[] = { 0, 1, 2, 3 };

template<>
double FE_P<3,1>::shape_value(const unsigned int i, const vec3 &p) const
{
    ASSERT(i <= 3, "Index of shape function is out of range.");
    switch (i)
    {
    case 0:
        return 1-p(0)-p(1)-p(2);
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
vec3 FE_P<3,1>::shape_grad(const unsigned int i, const vec3 &p) const
{
    ASSERT(i <= 3, "Index of shape function is out of range.");
    switch (i)
    {
    case 0:
        return "-1; -1; -1";
        break;
    case 1:
        return "1; 0; 0";
        break;
    case 2:
        return "0; 1; 0";
        break;
    case 3:
        return "0; 0; 1";
        break;
    }
}



#endif /* FE_P_HH_ */
