/*
 * fe_p.cc
 *
 *  Created on: Sep 3, 2012
 *      Author: jb
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



