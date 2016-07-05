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
 * @file    ref_element.cc
 * @brief   Class RefElement defines numbering of vertices, sides, calculation of normal vectors etc.
 * @author  Jan Stebel
 */

#include "system/global_defs.h"
#include "system/system.hh"
#include "mesh/ref_element.hh"



using namespace arma;
using namespace std;

template<unsigned int n>
std::vector< std::vector<unsigned int> > _array_to_vec( const unsigned int array[][n], unsigned int m) {
    std::vector< std::vector<unsigned int> > vec(m);
    for(unsigned int i=0; i<m; i++)
        for(unsigned int j=0;j<n; j++)
            vec[i].push_back(array[i][j]);
    return vec;
}

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


template<> const unsigned int RefElement<1>::side_nodes[][1] = {
		{ 0 },
		{ 1 }
};

template<> const unsigned int RefElement<2>::side_nodes[][2] = {
        { 0, 1},
        { 0, 2},
        { 1, 2}
};

template<> const unsigned int RefElement<3>::side_nodes[][3] = {
        {0,1,2},
        {0,1,3},
        {0,2,3},
        {1,2,3}
};



template<> const unsigned int RefElement<3>::side_lines[][3] = {
        {0,1,2},
        {0,3,4},
        {1,3,5},
        {2,4,5}
};



template<> const unsigned int RefElement<3>::line_nodes[][2] = {
        {0,1},
        {0,2},
        {1,2},
        {0,3},
        {1,3},
        {2,3}
};


template<> const std::vector< std::vector< std::vector<unsigned int> > > RefElement<1>::nodes_of_subelements = {
        _array_to_vec(side_nodes, n_sides),
        { {0,1} }
};

template<> const std::vector< std::vector< std::vector<unsigned int> > > RefElement<2>::nodes_of_subelements = {
        { {0}, {1}, {2} },
        _array_to_vec(side_nodes, n_sides),
        { {0,1,2} }
};

template<> const std::vector< std::vector< std::vector<unsigned int> > > RefElement<3>::nodes_of_subelements = {
        { {0}, {1}, {2}, {3} },
        _array_to_vec(line_nodes, n_lines),
        _array_to_vec(side_nodes, n_sides),
        { {0,1,2,3} }
};


template<unsigned int dim>
vec::fixed<dim> RefElement<dim>::node_coords(unsigned int nid)
{
	OLD_ASSERT(nid < n_nodes, "Vertex number is out of range!");

	vec::fixed<dim> p;
	p.zeros();

	if (nid > 0)
		p(nid-1) = 1;

	return p;
}


template<unsigned int dim>
vec::fixed<dim+1> RefElement<dim>::node_barycentric_coords(unsigned int nid)
{
	OLD_ASSERT(nid < n_nodes, "Vertex number is out of range!");

	vec::fixed<dim+1> p;
	p.zeros();

	if (nid == 0)
		p(dim) = 1;
	else
		p(nid-1) = 1;

	return p;
}


template<>
vec::fixed<1> RefElement<1>::normal_vector(unsigned int sid)
{
	OLD_ASSERT(sid < n_sides, "Side number is out of range!");

	return node_coords(sid) - node_coords(1-sid);
}

template<>
vec::fixed<2> RefElement<2>::normal_vector(unsigned int sid)
{
	OLD_ASSERT(sid < n_sides, "Side number is out of range!");
	vec::fixed<2> barycenter, bar_side, n, t;

	// tangent vector along line
	t = node_coords(side_nodes[sid][1]) - node_coords(side_nodes[sid][0]);
	// barycenter coordinates
	barycenter.fill(1./3);
	// vector from barycenter to the side
	bar_side = node_coords(side_nodes[sid][0]) - barycenter;
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
	OLD_ASSERT(sid < n_sides, "Side number is out of range!");
	vec::fixed<3> barycenter, bar_side, n, t1, t2;

	// tangent vectors of side
	t1 = node_coords(side_nodes[sid][1]) - node_coords(side_nodes[sid][0]);
	t2 = node_coords(side_nodes[sid][2]) - node_coords(side_nodes[sid][0]);
	// baryucenter coordinates
	barycenter.fill(0.25);
	// vector from barycenter to the side
	bar_side = node_coords(side_nodes[sid][0]) - barycenter;
	// normal vector (modulo sign)
	n = cross(t1,t2);
	n /= norm(n,2);
	// check sign of normal vector
	if (dot(n,bar_side) < 0) n = -n;

	return n;
}



template<unsigned int dim>
auto RefElement<dim>::barycentric_on_face(const BaryPoint &barycentric, unsigned int i_face) -> FaceBaryPoint
{
    ASSERT_EQ_DBG(barycentric.n_rows, dim+1);
    FaceBaryPoint face_barycentric;
    for(unsigned int i=0; i < dim; i++) {
        unsigned int i_sub_node = (i+1)%dim;
        unsigned int i_bary = (dim + side_nodes[i_face][i_sub_node])%(dim+1);
        face_barycentric[i] = barycentric[ i_bary ];
    }
    return face_barycentric;
}


template<unsigned int dim>
auto RefElement<dim>::barycentric_from_face(const FaceBaryPoint &face_barycentric, unsigned int i_face) -> BaryPoint
{
    ASSERT_EQ_DBG(face_barycentric.n_rows, dim);
    BaryPoint barycentric;
    barycentric.zeros();
    for(unsigned int i_sub_coord=0; i_sub_coord<dim; i_sub_coord++) {
        unsigned int i_sub_node = (i_sub_coord+1)%dim;
        barycentric+=face_barycentric(i_sub_coord)*node_barycentric_coords( side_nodes[i_face][i_sub_node]);
    }
    return barycentric;
}

template<>
auto RefElement<0>::clip(const BaryPoint &barycentric) -> BaryPoint
{
    return barycentric;
}

template<unsigned int dim>
auto  RefElement<dim>::make_bary_unit_vec()->BarycentricUnitVec
{
    std::vector<arma::vec::fixed<dim+1> > bary_unit_vec(dim+1, arma::zeros(dim+1));
    for(unsigned int i=0; i<dim; i++) {
        bary_unit_vec[i][i] = 1.0;
        bary_unit_vec[i][dim] = -1.0;
        bary_unit_vec[dim][i] = -1.0 / dim;
    }
    bary_unit_vec[dim][dim] = 1.0;
    return bary_unit_vec;
}


template<unsigned int dim>
auto RefElement<dim>::clip(const BaryPoint &barycentric) -> BaryPoint
{
    static BarycentricUnitVec bary_unit_vec = make_bary_unit_vec();
    ASSERT_EQ_DBG(barycentric.n_rows, dim+1);
    for(unsigned int i_bary=0; i_bary < dim +1; i_bary ++) {
        if (barycentric[i_bary] < 0.0) {
            // index of barycentric coord that is constant on the face i_side
            // as we use barycentric coords starting with local coordinates:
            // TODO: rather work only with local coords and/or with canonical barycentric coords
            unsigned int i_side = (2*dim - i_bary)%(dim +1);
            // project to face
            arma::vec projection_to_face(dim+1);
            //barycentric.print(cout, "input");
            //cout << "is: " << i_side << endl;
            //cout << "ibary: " << i_bary << endl;
            //bary_unit_vec[i_bary].print(cout, "normal");
            //barycentric.subvec(0, dim-1).print(cout, "bary sub");
            projection_to_face = barycentric - barycentric[i_bary]*bary_unit_vec[i_bary];
            //projection_to_face(dim) = 1.0 - arma::sum(projection_to_face.subvec(0, dim-1));
            //projection_to_face.print(cout, "projection");
            auto bary_on_face = barycentric_on_face(projection_to_face, i_side);
            //bary_on_face.print(cout, "b on f");
            auto sub_clip = RefElement<dim-1>::clip(bary_on_face);
            //sub_clip.print(cout, "sub clip");
            return barycentric_from_face(sub_clip, i_side);
        }
    }
    return barycentric;

}



template<unsigned int dim>
auto RefElement<dim>::centers_of_subelements(unsigned int sub_dim)->CentersList
{
    static std::vector< std::vector<LocalPoint> > list;
    if (list.size() == 0) {
        list.resize(dim+1);
        for(unsigned int sdim=0; sdim < dim+1; sdim++) {
            // Temporary solution until we have unified interface to
            // the internal indexing.
            // We use the fact that numbering of subelements goes as numbering of
            // k combinations over nodes.
            std::vector<unsigned int> subel_comb(sdim+2);
            for(auto &sub_el_nodes : nodes_of_subelements[sdim]) {
                ASSERT_EQ_DBG(sub_el_nodes.size(), sdim+1);
                LocalPoint center = arma::zeros(dim);
                for( unsigned int i_node : sub_el_nodes)
                    center+=node_coords( i_node );
                center/=(sdim+1);
                list[sdim].push_back(center);
            }
        }
    }

    ASSERT_LE_DBG(sub_dim, dim);
    return list[sub_dim];
}


template<>
double RefElement<1>::side_measure(unsigned int sid)
{
	OLD_ASSERT(sid < n_sides, "Side number is out of range!");

	return 1;
}


template<>
double RefElement<2>::side_measure(unsigned int sid)
{
	OLD_ASSERT(sid < n_sides, "Side number is out of range!");

	return norm(node_coords(side_nodes[sid][1]) - node_coords(side_nodes[sid][0]),2);
}


template<>
double RefElement<3>::side_measure(unsigned int sid)
{
	OLD_ASSERT(sid < n_sides, "Side number is out of range!");

	return 0.5*norm(cross(node_coords(side_nodes[sid][1]) - node_coords(side_nodes[sid][0]),
			node_coords(side_nodes[sid][2]) - node_coords(side_nodes[sid][0])),2);
}



template <>
unsigned int RefElement<3>::line_between_faces(unsigned int f1, unsigned int f2) {
    unsigned int i,j;
    i=j=0;
    while (side_lines[f1][i] != side_lines[f2][j])
        if (side_lines[f1][i] < side_lines[f2][j]) i++;
        else j++;
    return side_lines[f1][i];
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


template class RefElement<1>;
template class RefElement<2>;
template class RefElement<3>;




