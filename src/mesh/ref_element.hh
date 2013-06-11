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
 *
 *
 * TODO:
 *
 * - design interface in such a way, that we can change numbering
 * - design numbering and orientations on ref element that is consistent (orientation and numbering od 2d el. match sides of 3d),
 *   and possibly allows permute vertices of elements so that sides sharing an edge match numbering and orientation (we dos'nt need permute faces)
 *
 *   Proposal(prefers combinatoric order) :
 *   1D  - orientation V0 -> V1
 *
 *   2D  - edges: E0: V0 -> V1,
 *                E1: V0 -> V2
 *                E2: V1 -> V2
 *         This maximize number of edge orientations achievable by edge permutations
 *         edge numbering   edge orientation( in original numbering)
 *         0 1 2            + + +
 *         0 2 1            - + +
 *         1 0 2            + + -
 *         1 2 0            - - +
 *         2 0 1            + - -
 *         2 1 0            - - -
 *
 *                   vertices   edges       normal (out = +)
 *   3D - sides: S0: 0 1 2      E0 E1 E3    -
 *               S1: 0 1 3      E0 E2 E4    +
 *               S2: 0 2 3      E1 E2 E5    -
 *               S3: 1 2 3      E3 E4 E5    -
 *
 *        edges: E0: V0 -> V1  x direction
 *               E1: V0 -> V2  y direction
 *               E2: V0 -> V3  z direction
 *               E3: V1 -> V2
 *               E4: V1 -> V3
 *               E5: V2 -> V3
 *
 * - functions from DEAL.ii:
 *   bool is_inside_unit_cell( point )
 *   line_to_cell_vertices(line, vertex) vertex index on line to index on whole element
 *   face_to_cell_vertices(face, vertex, Orientation), Orientation should be some class describing permutation of shared face to element's face/side
 *   face_to_cel_lines
 *   standard_to_real_face_vertex(vertex, Orientation), maps vertex to permuted face
 *   real_to_standard_face_vertex - inverse
 *   ... same for line; we should need also something like standard_to_real_line_vertex; Deal dosn;t support Orientation changes for 2D element faces
 *   Point unit_cell_vertex(vertex) - coordinates
 *   project_to_unit_cell
 *   distance_to_unit_cell
 *   d_linear_shape_function
 *
 * - can not change numbering of element sides due to DarcyFlow, which use hardwired side numbering in construction of basis functions
 * - any change of side numbering requires also change in flow/old_bcd.cc
 *
 *
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
 * side id  node ids           side id  node ids           side id  node ids   normal
 * 0        1                  0        1,2                0        1,2,3      OUT
 * 1        0                  1        0,2                1        0,2,3      IN
 *                             2        0,1                2        0,1,3      OUT
 *                                                         3        0,1,2      IN
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
	 * Return index of 1D line, shared by two faces @p f1 and @p f2 of the reference tetrahedron.
	 * Implemented only for @p dim == 3.
	 */
	static unsigned int line_between_faces(unsigned int f1, unsigned int f2);


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

	/// Number of lines on boundary of one side.
	static const unsigned int n_lines_per_side = ( dim == 3 ? 3 : 0);

	/// Number of lines, i.e. @p object of dimension @p dim-2 on the boundary of the reference element.
	static const unsigned int n_lines = ( dim == 3 ? 6 : 0);

	/**
	 * Node numbers for each side.
	 */
	static const unsigned int side_nodes[n_sides][n_nodes_per_side];

	/**
	 * Indices of 1D lines of the 2D sides of an tetrahedron. Nonempty only for @p dim==3.
	 */
	static const unsigned int side_lines[n_sides][n_lines_per_side];

	/**
	 * Nodes of 1D lines of the tetrahedron.
	 */
    static const unsigned int line_nodes[n_lines][2];

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
