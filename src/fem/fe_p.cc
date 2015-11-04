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
 * @file    fe_p.cc
 * @brief   
 */

// !! implementation of specializations has to be i *.cc file to avoid multiple definition error during linking
#include "fe_p.hh"
#include "mesh/ref_element.hh"

/****** Template specializations ******/

/*** 1D finite elements ***/

// P0 constant element
template<>
DofDistribution<0,1>::DofDistribution()
{
    number_of_dofs = 1;

    number_of_single_dofs[1] = 1;

    unit_support_points.push_back(arma::zeros<arma::vec>(1));
}

// P1 linear element
template<>
DofDistribution<1,1>::DofDistribution()
{
    number_of_dofs = 2;

    number_of_single_dofs[0] = 2;

    unit_support_points.push_back(RefElement<1>::node_coords(0));
    unit_support_points.push_back(RefElement<1>::node_coords(1));
}

// P2 quadratic element
template<>
DofDistribution<2,1>::DofDistribution()
{
    number_of_dofs = 3;

    number_of_single_dofs[0] = 2;
    number_of_single_dofs[1] = 1;

    unit_support_points.push_back(RefElement<1>::node_coords(0));
    unit_support_points.push_back(RefElement<1>::node_coords(1));
    unit_support_points.push_back(
    		(RefElement<1>::node_coords(0)+RefElement<1>::node_coords(1))*0.5);
}

// P3 cubic element
template<>
DofDistribution<3,1>::DofDistribution()
{
    number_of_dofs = 4;

    number_of_single_dofs[0] = 2;
    number_of_pairs[1] = 1;

    unit_support_points.push_back(RefElement<1>::node_coords(0));
    unit_support_points.push_back(RefElement<1>::node_coords(1));
    unit_support_points.push_back(
    		(RefElement<1>::node_coords(0)*2+RefElement<1>::node_coords(1))/3);
    unit_support_points.push_back(
    		(RefElement<1>::node_coords(0)+RefElement<1>::node_coords(1)*2)/3);
}


/*** 2D finite elements ***/

// P0 constant element
template<>
DofDistribution<0,2>::DofDistribution()
{
    number_of_dofs = 1;

    number_of_single_dofs[2] = 1;

    unit_support_points.push_back(RefElement<2>::node_coords(0));
}


// P1 linear element
template<>
DofDistribution<1,2>::DofDistribution()
{
    number_of_dofs = 3;

    number_of_single_dofs[0] = 3;

    unit_support_points.push_back(RefElement<2>::node_coords(0));
    unit_support_points.push_back(RefElement<2>::node_coords(1));
    unit_support_points.push_back(RefElement<2>::node_coords(2));
}

// P2 quadratic element
template<>
DofDistribution<2,2>::DofDistribution()
{
    number_of_dofs = 6;

    number_of_single_dofs[0] = 3;
    number_of_single_dofs[1] = 3;

    unit_support_points.push_back(RefElement<2>::node_coords(0));
    unit_support_points.push_back(RefElement<2>::node_coords(1));
    unit_support_points.push_back(RefElement<2>::node_coords(2));
    unit_support_points.push_back(
    		(RefElement<2>::node_coords(0)+RefElement<2>::node_coords(1))*0.5);
    unit_support_points.push_back(
    		(RefElement<2>::node_coords(0)+RefElement<2>::node_coords(2))*0.5);
    unit_support_points.push_back(
    		(RefElement<2>::node_coords(1)+RefElement<2>::node_coords(2))*0.5);
}

// P3 cubic element
template<>
DofDistribution<3,2>::DofDistribution()
{
    number_of_dofs = 10;

    number_of_single_dofs[0] = 3;
    number_of_pairs[1] = 3;
    number_of_single_dofs[2] = 1;

    unit_support_points.push_back(RefElement<2>::node_coords(0));
    unit_support_points.push_back(RefElement<2>::node_coords(1));
    unit_support_points.push_back(RefElement<2>::node_coords(2));
    unit_support_points.push_back(
    		(RefElement<2>::node_coords(0)*2+RefElement<2>::node_coords(1))/3);
    unit_support_points.push_back(
    		(RefElement<2>::node_coords(0)+RefElement<2>::node_coords(1)*2)/3);
    unit_support_points.push_back(
    		(RefElement<2>::node_coords(0)*2+RefElement<2>::node_coords(2))/3);
    unit_support_points.push_back(
    		(RefElement<2>::node_coords(0)+RefElement<2>::node_coords(2)*2)/3);
    unit_support_points.push_back(
    		(RefElement<2>::node_coords(1)*2+RefElement<2>::node_coords(2))/3);
    unit_support_points.push_back(
    		(RefElement<2>::node_coords(1)+RefElement<2>::node_coords(2)*2)/3);
    unit_support_points.push_back(
    		(RefElement<2>::node_coords(0)+RefElement<2>::node_coords(1)+RefElement<2>::node_coords(2))/3);
}



/*** 3D finite elements ***/

// P0 constant element
template<>
DofDistribution<0,3>::DofDistribution()
{
    number_of_dofs = 1;

    number_of_single_dofs[3] = 1;

    unit_support_points.push_back(RefElement<3>::node_coords(0));
}


// P1 linear element
template<>
DofDistribution<1,3>::DofDistribution()
{
    number_of_dofs = 4;

    number_of_single_dofs[0] = 4;

    unit_support_points.push_back(RefElement<3>::node_coords(0));
    unit_support_points.push_back(RefElement<3>::node_coords(1));
    unit_support_points.push_back(RefElement<3>::node_coords(2));
    unit_support_points.push_back(RefElement<3>::node_coords(3));
}

// P2 quadratic element
template<>
DofDistribution<2,3>::DofDistribution()
{
    number_of_dofs = 10;

    number_of_single_dofs[0] = 4;
    number_of_single_dofs[1] = 6;

    unit_support_points.push_back(RefElement<3>::node_coords(0));
    unit_support_points.push_back(RefElement<3>::node_coords(1));
    unit_support_points.push_back(RefElement<3>::node_coords(2));
    unit_support_points.push_back(RefElement<3>::node_coords(3));
    unit_support_points.push_back(
    		(RefElement<3>::node_coords(0)+RefElement<3>::node_coords(1))*0.5);
    unit_support_points.push_back(
    		(RefElement<3>::node_coords(0)+RefElement<3>::node_coords(2))*0.5);
    unit_support_points.push_back(
    		(RefElement<3>::node_coords(0)+RefElement<3>::node_coords(3))*0.5);
    unit_support_points.push_back(
    		(RefElement<3>::node_coords(1)+RefElement<3>::node_coords(2))*0.5);
    unit_support_points.push_back(
    		(RefElement<3>::node_coords(1)+RefElement<3>::node_coords(3))*0.5);
    unit_support_points.push_back(
    		(RefElement<3>::node_coords(2)+RefElement<3>::node_coords(3))*0.5);
}

// P3 cubic element
template<>
DofDistribution<3,3>::DofDistribution()
{
    number_of_dofs = 20;

    number_of_single_dofs[0] = 4;
    number_of_pairs[1] = 6;
    number_of_single_dofs[2] = 4;

    unit_support_points.push_back(RefElement<3>::node_coords(0));
    unit_support_points.push_back(RefElement<3>::node_coords(1));
    unit_support_points.push_back(RefElement<3>::node_coords(2));
    unit_support_points.push_back(RefElement<3>::node_coords(3));
    unit_support_points.push_back(
    		(RefElement<3>::node_coords(0)*2+RefElement<3>::node_coords(1))/3);
    unit_support_points.push_back(
    		(RefElement<3>::node_coords(0)+RefElement<3>::node_coords(1)*2)/3);
    unit_support_points.push_back(
    		(RefElement<3>::node_coords(0)*2+RefElement<3>::node_coords(2))/3);
    unit_support_points.push_back(
    		(RefElement<3>::node_coords(0)+RefElement<3>::node_coords(2)*2)/3);
    unit_support_points.push_back(
    		(RefElement<3>::node_coords(0)*2+RefElement<3>::node_coords(3))/3);
    unit_support_points.push_back(
    		(RefElement<3>::node_coords(0)+RefElement<3>::node_coords(3)*2)/3);
    unit_support_points.push_back(
    		(RefElement<3>::node_coords(1)+RefElement<3>::node_coords(2)*2)/3);
    unit_support_points.push_back(
    		(RefElement<3>::node_coords(1)*2+RefElement<3>::node_coords(2))/3);
    unit_support_points.push_back(
    		(RefElement<3>::node_coords(1)+RefElement<3>::node_coords(3)*2)/3);
    unit_support_points.push_back(
    		(RefElement<3>::node_coords(1)*2+RefElement<3>::node_coords(3))/3);
    unit_support_points.push_back(
    		(RefElement<3>::node_coords(2)+RefElement<3>::node_coords(3)*2)/3);
    unit_support_points.push_back(
    		(RefElement<3>::node_coords(2)*2+RefElement<3>::node_coords(3))/3);
    unit_support_points.push_back(
    		(RefElement<3>::node_coords(0)+RefElement<3>::node_coords(1)+RefElement<3>::node_coords(2))/3);
    unit_support_points.push_back(
    		(RefElement<3>::node_coords(0)+RefElement<3>::node_coords(1)+RefElement<3>::node_coords(3))/3);
    unit_support_points.push_back(
    		(RefElement<3>::node_coords(0)+RefElement<3>::node_coords(2)+RefElement<3>::node_coords(3))/3);
    unit_support_points.push_back(
    		(RefElement<3>::node_coords(1)+RefElement<3>::node_coords(2)+RefElement<3>::node_coords(3))/3);
}



