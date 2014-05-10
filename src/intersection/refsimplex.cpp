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
#include "refsimplex.h"
#include <array>

using namespace arma;
using namespace std;

namespace computeintersection{

template<> const unsigned int RefSimplex<1>::side_permutations[][n_nodes_per_side] = { { 0 } };

template<> const unsigned int RefSimplex<2>::side_permutations[][n_nodes_per_side] = { { 0, 1 }, { 1, 0 } };

template<> const unsigned int RefSimplex<3>::side_permutations[][n_nodes_per_side] = {
		{ 0, 1, 2 },
		{ 0, 2, 1 },
		{ 1, 0, 2 },
		{ 1, 2, 0 },
		{ 2, 0, 1 },
		{ 2, 1, 0 }
};


template<> const unsigned int RefSimplex<1>::side_nodes[][1] = {
		{ 0 },
		{ 1 }
};

template<> const unsigned int RefSimplex<2>::side_nodes[][2] = {
        { 0, 1},
        { 0, 2},
        { 1, 2}
};

template<> const unsigned int RefSimplex<3>::side_nodes[][3] = {
        {0,1,2},
        {0,1,3},
        {0,2,3},
        {1,2,3}
};



template<> const unsigned int RefSimplex<3>::side_lines[][3] = {
		{0,1,2},
		{0,3,4},
		{1,3,5},
		{2,4,5}
};

//template<unsigned int dim>
//const unsigned int RefElement<dim>::side_lines[][0] = {{}};


/*
 * Indexes of nodes for each line
 * */
template<> const unsigned int RefSimplex<3>::line_nodes[][2] = {
        {0,1},
        {0,2},
        {1,2},
        {0,3},
        {1,3},
        {2,3}
};

template<> const unsigned int RefSimplex<2>::line_nodes[][2] = {
		{0,1},
		{0,2},
		{1,2}
};

/**
 * Indexes of sides for each line - with right orientation
 */

template<> const unsigned int RefSimplex<3>::line_sides[][2] = {
		{0,1},
		{2,0},
		{0,3},
		{1,2},
		{3,1},
		{2,3}
};


//template<unsigned int dim>
//const unsigned int RefElement<dim>::line_nodes[][0] = {};


//template<> static std::array< arma::vec::fixed<dim+1>, subdim+1 > bary_coords(unsigned int sid){

//};
template<unsigned int dim>
template<unsigned int subdim>
RefSimplex<subdim> RefSimplex<dim>::SubRefSimplex(){

};

/*template<unsigned int dim>
template<unsigned int subdim>
std::array< arma::vec::fixed<dim+1>,subdim+1> RefSimplex<dim>::bary_coords(){

	//ASSERT(subdim < dim, "Sub-dimension is bigger than dimension!");
	xprintf(Msg, "barycoods \n");

	std::array<arma::vec::fixed<dim+1>,subdim+1> bary_c;
	for(unsigned int i = 0; i < dim; i++){
		bary_c[i] = RefSimplex<dim>::node_coords(RefSimplex<dim>::side_nodes[sid][i]);
	}

	return bary_c;
};*/


template<unsigned int dim>
vec::fixed<dim+1> RefSimplex<dim>::node_coords(unsigned int nid)
{
	ASSERT(nid < n_nodes, "Vertex number is out of range!");

	vec::fixed<dim+1> p;
	p.zeros();

	p(nid) = 1;

	return p;
}


template<unsigned int dim>
vec::fixed<dim> RefSimplex<dim>::normal_vector(unsigned int sid)
{
	ASSERT(sid < n_sides, "Side number is out of range!");
	vec::fixed<dim> p;
	unsigned int new_sid = sid;

	if (dim==1)
		new_sid = (sid+1)%2;
	else if (dim==2)
		new_sid = (sid+2)%3;

	if (new_sid == 0)
		p.fill(1./sqrt(dim));
	else
	{
		p.zeros();
		p(new_sid-1) = -1;
	}

	return p;
}



template <>
unsigned int RefSimplex<3>::line_between_faces(unsigned int f1, unsigned int f2) {
    unsigned int i,j;
    i=j=0;
    while (side_lines[f1][i] != side_lines[f2][j])
        if (side_lines[f1][i] < side_lines[f2][j]) i++;
        else j++;
    return side_lines[f1][i];
}



template<unsigned int dim>
unsigned int RefSimplex<dim>::permutation_index(unsigned int p[n_nodes_per_side])
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


template class RefSimplex<1>;
template class RefSimplex<2>;
template class RefSimplex<3>;

} // END namespace


