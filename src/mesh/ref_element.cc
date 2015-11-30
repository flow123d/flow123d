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
 * @brief Class RefElement defines numbering of vertices, sides, calculation of normal vectors etc.
 * @author Jan Stebel
 */

#include "system/global_defs.h"
#include "system/system.hh"
#include "mesh/ref_element.hh"

using namespace arma;
using namespace std;

    
    
template<> const IdxVector<2> RefElement<1>::line_nodes_[] = {
        {0,1}
};

template<> const IdxVector<2> RefElement<2>::line_nodes_[] = {
        {0,1},
        {0,2},
        {1,2}
};

template<> const IdxVector<2> RefElement<3>::line_nodes_[] = {
        {0,1},
        {0,2},
        {1,2},
        {0,3},
        {1,3},
        {2,3}
};


template<> const IdxVector<1> RefElement<1>::node_lines_[] = {
     {0},
     {0}
};

template<> const IdxVector<2> RefElement<2>::node_lines_[] = {
     {1,0},
     {0,2},
     {2,1}
};

//TODO: what should be the order ??
template<> const IdxVector<3> RefElement<3>::node_lines_[] = {
     {0,1,3},
     {0,2,4},
     {1,2,5},
     {3,4,5},
};


//TODO: what should be the order ??
template<> const IdxVector<3> RefElement<3>::side_nodes_[] = {
        { 0, 1, 2 },
        { 0, 1, 3 },
        { 0, 2, 3 },
        { 1, 2, 3 }
};

//TODO: what should be the order ??
template<> const IdxVector<3> RefElement<3>::node_sides_[] = {
        { 0, 1, 2 },
        { 0, 1, 3 },
        { 0, 2, 3 },
        { 1, 2, 3 }
};

template<> const IdxVector<2> RefElement<3>::line_sides_[] = {
     {0,1},
     {2,0},
     {0,3},
     {1,2},
     {3,1},
     {2,3}
};


template<> const IdxVector<3> RefElement<3>::side_lines_[] = {
        {0,1,2},
        {0,3,4},
        {1,3,5},
        {2,4,5}
};



template<> const unsigned int RefElement<1>::side_permutations[][n_nodes_per_side] = { { 0 } };

template<> const unsigned int RefElement<2>::side_permutations[][n_nodes_per_side] = { { 0, 1 }, { 1, 0 } };

template<> const unsigned int RefElement<3>::side_permutations[][n_nodes_per_side] = {
		{ 0, 1, 2 },
		{ 0, 2, 1 },
		{ 1, 0, 2 },
		{ 1, 2, 0 },
		{ 2, 0, 1 },
		{ 2, 1, 0 }
};


// template<> const unsigned int RefElement<1>::side_nodes[][1] = {
// 		{ 0 },
// 		{ 1 }
// };
// 
// template<> const unsigned int RefElement<2>::side_nodes[][2] = {
// 		{ 0, 1 },
// 		{ 0, 2 },
// 		{ 1, 2 }
// };
// 
// template<> const unsigned int RefElement<3>::side_nodes[][3] = {
// 		{ 0, 1, 2 },
// 		{ 0, 1, 3 },
// 		{ 0, 2, 3 },
// 		{ 1, 2, 3 }
// };
// 
// 
// 
// template<> const unsigned int RefElement<3>::side_lines[][3] = {
//         {0,1,2},
//         {0,3,4},
//         {1,3,5},
//         {2,4,5}
// };
// 
// 
// template<> const unsigned int RefElement<1>::line_nodes[][2] = {
//         {0,1}
// };
// 
// template<> const unsigned int RefElement<2>::line_nodes[][2] = {
//         {0,1},
//         {0,2},
//         {1,2}
// };
// 
// template<> const unsigned int RefElement<3>::line_nodes[][2] = {
//         {0,1},
//         {0,2},
//         {1,2},
//         {0,3},
//         {1,3},
//         {2,3}
// };
// 
// 
// /**
//  * Indexes of sides for each line - with right orientation
//  */
// 
// template<> const unsigned int RefElement<3>::line_sides[][2] = {
//      {0,1},
//      {2,0},
//      {0,3},
//      {1,2},
//      {3,1},
//      {2,3}
// };
// 
// /**
//  * Indexes of sides for each line - for Simplex<2>, with right orientation
//  */
// template<> const unsigned int RefElement<2>::line_sides[][2] = {
//      {1,0},
//      {0,2},
//      {2,1}
// };


template<unsigned int dim>
vec::fixed<dim> RefElement<dim>::node_coords(unsigned int nid)
{
	ASSERT(nid < n_nodes, "Vertex number is out of range!");

	vec::fixed<dim> p;
	p.zeros();

    // these are real coordinates in x,y,z as we want
    if (nid > 0)
        p(nid-1) = 1;

	return p;
}


template<unsigned int dim>
vec::fixed<dim+1> RefElement<dim>::node_barycentric_coords(unsigned int nid)
{
    ASSERT(nid < n_nodes, "Vertex number is out of range!");

    vec::fixed<dim+1> p;
    p.zeros();

// this is by VF    
    p(nid) = 1;

    // this is what we want
//     if (nid == 0)
//         p(dim) = 1;
//     else
//         p(nid-1) = 1;
    
    return p;
}


template<unsigned int dim>
inline unsigned int RefElement<dim>::oposite_node(unsigned int sid)
{
    return n_sides - sid - 1;
}

template<>
unsigned int RefElement<3>::normal_orientation(unsigned int sid)
{
    ASSERT(sid < n_sides, "Side number is out of range!");

    return (sid+1) % 2;
}

template<unsigned int dim>
unsigned int RefElement<dim>::normal_orientation(unsigned int sid)
{
    ASSERT(sid < n_sides, "Side number is out of range!");

    return sid % 2;
}


template<>
vec::fixed<1> RefElement<1>::normal_vector(unsigned int sid)
{
    ASSERT(sid < n_sides, "Side number is out of range!");

    return node_coords(sid) - node_coords(1-sid);
}

template<>
vec::fixed<2> RefElement<2>::normal_vector(unsigned int sid)
{
    ASSERT(sid < n_sides, "Side number is out of range!");
    vec::fixed<2> barycenter, bar_side, n, t;

    // tangent vector along line
    t = node_coords(line_nodes_[sid][1]) - node_coords(line_nodes_[sid][0]);
    // barycenter coordinates
    barycenter.fill(1./3);
    // vector from barycenter to the side
    bar_side = node_coords(line_nodes_[sid][0]) - barycenter;
    // normal vector to side (modulo sign)
    n(0) = -t(1);
    n(1) = t(0);
    n /= norm(n,2);
    // check sign of normal vector
    if (dot(n,bar_side) < 0) n *= -1;

    return n;
}

template<>
vec::fixed<3> RefElement<3>::normal_vector(unsigned int sid)
{
    ASSERT(sid < n_sides, "Side number is out of range!");
    vec::fixed<3> barycenter, bar_side, n, t1, t2;

    // tangent vectors of side
    t1 = node_coords(side_nodes_[sid][1]) - node_coords(side_nodes_[sid][0]);
    t2 = node_coords(side_nodes_[sid][2]) - node_coords(side_nodes_[sid][0]);
    // baryucenter coordinates
    barycenter.fill(0.25);
    // vector from barycenter to the side
    bar_side = node_coords(side_nodes_[sid][0]) - barycenter;
    // normal vector (modulo sign)
    n = cross(t1,t2);
    n /= norm(n,2);
    // check sign of normal vector
    if (dot(n,bar_side) < 0) n = -n;

    return n;
}


template<>
double RefElement<1>::side_measure(unsigned int sid)
{
    ASSERT(sid < n_sides, "Side number is out of range!");

    return 1;
}


template<>
double RefElement<2>::side_measure(unsigned int sid)
{
    ASSERT(sid < n_sides, "Side number is out of range!");

    return norm(node_coords(line_nodes_[sid][1]) - node_coords(line_nodes_[sid][0]),2);
}


template<>
double RefElement<3>::side_measure(unsigned int sid)
{
    ASSERT(sid < n_sides, "Side number is out of range!");

    return 0.5*norm(cross(node_coords(side_nodes_[sid][1]) - node_coords(side_nodes_[sid][0]),
            node_coords(side_nodes_[sid][2]) - node_coords(side_nodes_[sid][0])),2);
}

template <>
unsigned int RefElement<3>::line_between_faces(unsigned int f1, unsigned int f2) {
    unsigned int i,j;
    i=j=0;
    while (side_lines_[f1][i] != side_lines_[f2][j])
        if (side_lines_[f1][i] < side_lines_[f2][j]) i++;
        else j++;
    return side_lines_[f1][i];
}



template<unsigned int dim>
unsigned int RefElement<dim>::permutation_index(unsigned int p[n_nodes_per_side])
{
	unsigned int index;
	for (index = 0; index < n_side_permutations; index++)
		if (equal(p, p + n_nodes_per_side, side_permutations[index]))
			return index;

	xprintf(PrgErr, "Side permutation not found.\n");

	// The following line is present in order to suppress compilers warning
	// about missing return value.
	return 0;
}


/**
     * Basic line interpolation
     */
template<unsigned int dim>
arma::vec::fixed<dim+1> RefElement<dim>::line_barycentric_interpolation(
                                                       arma::vec::fixed<dim+1> first_coords, 
                                                       arma::vec::fixed<dim+1> second_coords, 
                                                       double first_theta, double second_theta, double theta){

    arma::vec::fixed<dim+1> bary_interpolated_coords;
    bary_interpolated_coords = ((theta - first_theta) * second_coords + (second_theta - theta) * first_coords)
                               /(second_theta - first_theta);
    return bary_interpolated_coords;
};

template class RefElement<1>;
template class RefElement<2>;
template class RefElement<3>;


