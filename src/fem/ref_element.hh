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

#ifndef REF_ELEMENT_HH_
#define REF_ELEMENT_HH_

#include <armadillo>

/*
 * Ordering of nodes and sides in reference elements
 * =================================================
 *
 * 1D element (line segment)   2D element (triangle)        3D element (tetrahedron)
 *
 *                                                                            y
 *                                                                          .
 *                                                                        ,/
 *                                                                       /
 *                                                                    2
 *                             y                                    ,/|`\
 *                             ^                                  ,/  |  `\
 *                             |                                ,/    '.   `\
 *                             2                              ,/       |     `\
 *                             |`\                          ,/         |       `\
 *                             |  `\                       0-----------'.--------1 --> x
 *                             |    `\                      `\.         |      ,/
 *                             |      `\                       `\.      |    ,/
 *                             |        `\                        `\.   '. ,/
 * 0----------1 --> x          0----------1 --> x                    `\. |/
 *                                                                      `3
 *                                                                         `\.
 *                                                                            ` z
 *
 * side id  node ids           side id  node ids           side id  node ids
 * 0        1                  0        1,2                0        1,2,3
 * 1        0                  1        0,2                1        0,2,3
 *                             2        0,1                2        0,1,3
 *                                                         3        0,1,2
 *
 *
 *
 *
 *
 *
 */
template<unsigned int dim>
class RefElement
{
public:

	/**
	 * Return barycentric coordinates of given node.
	 * @param nid Node number.
	 */
	static arma::vec::fixed<dim+1> node_coords(unsigned int nid);

	/**
	 * Compute normal vector to a given side.
	 * @param sid Side number.
	 */
	static arma::vec::fixed<dim> normal_vector(unsigned int sid);


	/**
	 * Number of sides.
	 */
	static const unsigned int n_sides = dim + 1;

	/**
	 * Number of vertices.
	 */
	static const unsigned int n_nodes = dim + 1;

	/**
	 * Number of nodes on one side.
	 */
	static const unsigned int n_nodes_per_side = dim;

	/**
	 * Node numbers for each side.
	 */
	static const unsigned int side_nodes[n_sides][n_nodes_per_side];

	/**
	 * Number of permutations of nodes on sides.
	 * dim   value
	 * -----------
	 * 1     1
	 * 2     2
	 * 3     6
	 */
	static const unsigned int n_side_permutations = (dim+1)*(2*dim*dim-5*dim+6)/6;

	/**
	 * Permutations of nodes on sides.
	 */
	static const unsigned int side_permutations[n_side_permutations][n_nodes_per_side];

	/**
	 * For a given permutation @p p of nodes finds its index within @p side_permutations.
	 * @param p Permutation of nodes.
	 */
	static unsigned int permutation_index(unsigned int p[n_nodes_per_side]);

};







#endif /* REF_ELEMENT_HH_ */
