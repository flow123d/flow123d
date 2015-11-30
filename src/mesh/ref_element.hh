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
 * TODO:  reconsider following whether it is actual...
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
#include "system/system.hh"

/*
 * Ordering of nodes and sides in reference elements
 * =================================================
 *
 * TODO we want the following (22.10.):
 * 
 * 1D element (line segment)   2D element (triangle)        3D element (tetrahedron)
 *
 *                                                                            z
 *                                                                          .
 *                                                                        ,/
 *                                                                       /
 *                                                                    3
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
 *                                                                      `2
 *                                                                         `\.
 *                                                                            `y
 *
 * side id  node ids           side id  node ids           side id  node ids   normal
 * 0        0                  0        0,1                0        0,1,2      IN
 * 1        1                  1        0,2                1        0,1,3      OUT
 *                             2        1,2                2        0,2,3      IN
 *                                                         3        1,2,3      OUT
 *
 *
 * nodes coordinates:                                                               
 * 0        [0]                0        [0,0]              0        [0,0,0]    
 * 1        [1]                1        [1,0]              1        [1,0,0]    
 *                             2        [0,1]              2        [0,1,0]       
 *                                                         3        [0,0,1]    
 * 
 * barycentric coordinates of nodes:
 * 0        [0,1]              0        [0,0,1]            0        [0,0,0,1]
 * 1        [1,0]              1        [1,0,0]            1        [1,0,0,0]
 *                             2        [0,1,0]            2        [0,1,0,0]
 *                                                         3        [0,0,1,0]
 */

/** Auxilliary class representing vector of indices (unsigned int).
 * @tparam Size is the fixed size of the vector.
 */
template<unsigned int Size>
class IdxVector{
    unsigned int data_[Size];   ///< Array with indices.
    
    public:
        /// Constructor taking in array of indices.
        IdxVector(std::array<unsigned int,Size> data_in);
        /// Constructor enabling creating object with initializer list {...}.
        IdxVector(std::initializer_list<unsigned int> data_in);
        /// Getter for index @p idx.
        unsigned int operator[](unsigned int idx) const;
};


template<unsigned int dim>
class RefElement
{
public:
        
	/**
	 * Return coordinates of given node.
     * @see the class documentation @p RefElement
	 * @param nid Node number.
     * NOTE: Implementation is dependent on current node and side numbering.
	 */
	static arma::vec::fixed<dim> node_coords(unsigned int nid);
    
    /**
     * Return barycentric coordinates of given node.
     * @see the class documentation @p RefElement
     * @param nid Node number.
     */
    static arma::vec::fixed<dim+1> node_barycentric_coords(unsigned int nid);
    
	/**
	 * Compute normal vector to a given side.
	 * @param sid Side number.
	 */
	static arma::vec::fixed<dim> normal_vector(unsigned int sid);

    /** Returns orientation of the normal of side @p sid. 0 -> OUT, 1 -> IN.
     * NOTE: Implementation is dependent on current node and side numbering.
     */
    static unsigned int normal_orientation(unsigned int sid);
    
	static double side_measure(unsigned int sid);
    
    /**
     * Returns index of the node that is oposite to side of given index @p sid.
     * Note: It is dependent on current node and side numbering.
     * @param sid Side number.
     * NOTE: Implementation is dependent on current node and side numbering.
     */
    static unsigned int oposite_node(unsigned int sid);

	/**
	 * Return index of 1D line, shared by two faces @p f1 and @p f2 of the reference tetrahedron.
	 * Implemented only for @p dim == 3.
	 */
	static unsigned int line_between_faces(unsigned int f1, unsigned int f2);

    
	static const unsigned int n_sides = dim + 1;            ///< Number of sides.
	static const unsigned int n_nodes = dim + 1;            ///< Number of nodes.
	static const unsigned int n_nodes_per_side = dim;       ///< Number of nodes on one side.
    static const unsigned int n_lines_per_node = dim;       ///< Number of lines with one common node.
    static const unsigned int n_nodes_per_line = 2;         ///< Number of nodes in one line.
    static const unsigned int n_sides_per_line = 2;         ///< Number of sides with one common line. @p dim == 3.
    
	/// Number of lines on boundary of one side.
	static const unsigned int n_lines_per_side = (unsigned int)((dim * (dim - 1)) / 2);//( dim == 3 ? 3 : 0);// Kombinační číslo dim nad dvěma

	/// Number of lines, i.e. @p object of dimension @p dim-2 on the boundary of the reference element.
	static const unsigned int n_lines = (unsigned int)((dim * (dim + 1)) / 2); //( dim == 3 ? 6 : dim == 2 ? 3 : dim == 1 ? 1 : 0); součet posloupnosti

    
// 	/**
// 	 * Node numbers for each side.
// 	 */
// 	static const unsigned int side_nodes[n_sides][n_nodes_per_side];
// 
// 	/**
// 	 * Indices of 1D lines of the 2D sides of an tetrahedron. Nonempty only for @p dim==3.
// 	 */
// 	static const unsigned int side_lines[n_sides][n_lines_per_side];
// 
// 	/**
// 	 * Nodes of 1D lines of the tetrahedron.
// 	 */
//     static const unsigned int line_nodes[n_lines][2];
//     
//     /**
//      * Indices of sides for each line. Nonempty only for @p dim==3 and @p dim==2.
//      */
//     static const unsigned int line_sides[n_lines][2];

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
    
    /**
     * @param sid - index of a sub-simplex in a simplex
     * return an array of barycentric coordinates on <dim> simplex from <subdim> simplex
     * for example: simplex<3> - ABCD and its subsubsimplex<1> AD (line index: 3)
     * AD has barycoords for A (1,0), for D (0,1), but A in ABCD is (1,0,0,0) and D is (0,0,0,1)
     * this method creates array ((1,0,0,0),(0,0,0,1))
     */
    template<unsigned int subdim> static std::array<arma::vec::fixed<dim+1>,subdim+1> bary_coords(unsigned int sid);

    /** Interpolate barycentric coords to a higher dimension of a simplex.
     * @param coord - barycentric coords of a point on a sub-simplex
     * @param sub_simplex_idx - id of sub-simplex on a simplex
     */
    template<unsigned int subdim> static arma::vec::fixed<dim+1> interpolate(arma::vec::fixed<subdim+1> coord, int sub_simplex_idx);

    /**
     * Basic line interpolation.
     */
    static arma::vec::fixed<dim+1> line_barycentric_interpolation(arma::vec::fixed<dim+1> first_coords, 
                                                                  arma::vec::fixed<dim+1> second_coords, 
                                                                  double first_theta, double second_theta, double theta);

    /**
     * This method serves as an interface to topology information of the reference element.
     * It returns indices of OutDim-dimensional object
     * of InDim-dimnesional object of given index
     * in dim-dimnesional reference element.
     * @tparam OutDim - output dimension (give me node-0, line-1, side-2), <= dim
     * @tparam InDim - input dimension (for node-0, line-1, side-2), <= dim
     * @return vector of indices represented by @p IdxVector object.
     * 
     * possible calls:
     *  dim    OutDim  InDim  return
     * 1,2,3   0       1      InDim+1   - give me indices of nodes of line of given index
     *   3     0       2      InDim+1   - give me indices of nodes of a side (triangle) of given index
     *   3     1       2      InDim+1   - give me indices of lines of side (triangle) of given index
     *                               
     * 1,2,3   1       0     dim-InDim  - give me indices of lines with common node of given index
     *   3     2       0     dim-InDim  - give me indices of sides (triangles) with common node of given index
     *   3     2       1     dim-InDim  - give me indices of sides (triangles) with common line of given index 
     * 
     */
    template<unsigned int OutDim, unsigned int InDim> 
    static const IdxVector< (InDim>OutDim ? InDim+1 : dim-InDim) > interact(unsigned int index);

private:
    static const IdxVector<n_nodes_per_line> line_nodes_[n_lines]; ///< For given line, returns its nodes indices.
    static const IdxVector<n_lines_per_node> node_lines_[n_nodes]; ///< For given node, returns lines indices.
    static const IdxVector<n_nodes_per_side> side_nodes_[n_sides]; ///< For given side, returns nodes indices. For @p dim == 3.
    static const IdxVector<n_nodes_per_side> node_sides_[n_nodes]; ///< For given node, returns sides indices. For @p dim == 3.
    static const IdxVector<n_sides_per_line> line_sides_[n_lines]; ///< For given line, returns sides indices. For @p dim == 3.
    static const IdxVector<n_lines_per_side> side_lines_[n_sides]; ///< For given side, returns lines indices. For @p dim == 3.
};




/************************* template implementation ****************************/

template<unsigned int dim>
template<unsigned int subdim> 
std::array<arma::vec::fixed<dim+1>,subdim+1> RefElement<dim>::bary_coords(unsigned int sid){
        ASSERT(subdim < dim, "Dimension mismatch!");
        std::array<arma::vec::fixed<dim+1>,subdim+1> bary_c;

        for(unsigned int i = 0; i < subdim+1; i++)
            bary_c[i] = RefElement<dim>::node_barycentric_coords(RefElement<dim>::interact<0,subdim>(sid)[i]);

        return bary_c;
};

template<unsigned int dim>
template<unsigned int subdim> 
arma::vec::fixed<dim+1> RefElement<dim>::interpolate(arma::vec::fixed<subdim+1> coord, int sub_simplex_idx){
    std::array<arma::vec::fixed<dim+1>, subdim+1> simplex_M_vertices = RefElement<dim>::bary_coords<subdim>(sub_simplex_idx);
    arma::vec::fixed<dim+1> sum;
    sum.zeros();
    for(int i=0; i<subdim+1; i++) sum += coord[i]*simplex_M_vertices[i];
    return sum;
};
    


template <unsigned int Size>
IdxVector<Size>::IdxVector(std::array<unsigned int,Size> data_in)
: data_(data_in){}

template <unsigned int Size>
IdxVector<Size>::IdxVector(std::initializer_list<unsigned int> data_in)
{
    ASSERT(data_in.size() == Size, "Incorrect data size.");
    std::copy(data_in.begin(), data_in.end(), data_);
}

template <unsigned int Size>
inline unsigned int IdxVector<Size>::operator[](unsigned int idx) const
{   ASSERT(idx < Size, "Index out of bounds.");
    return data_[idx]; }
    



/// This function is for "side_nodes" - for given side, give me nodes (0->0, 1->1).
template<> template<> inline const IdxVector<1> RefElement<1>::interact<0,0>(unsigned int i)
{   ASSERT(i < RefElement<1>::n_nodes, "Index out of bounds.");
    return IdxVector<1>({i});}

/// For line i {0}, give me indices of its nodes.
template<> template<> inline const IdxVector<2> RefElement<1>::interact<0,1>(unsigned int i)
{   ASSERT(i < RefElement<1>::n_lines, "Index out of bounds.");
    return line_nodes_[i];}

/// For line i {0,1,2}, give me indices of its nodes.
template<> template<> inline const IdxVector<2> RefElement<2>::interact<0,1>(unsigned int i)
{   ASSERT(i < RefElement<2>::n_lines, "Index out of bounds.");
    return line_nodes_[i];}

/// For line i {0,1,2,3,4,5}, give me indices of its nodes.
template<> template<> inline const IdxVector<2> RefElement<3>::interact<0,1>(unsigned int i)
{   ASSERT(i < RefElement<3>::n_lines, "Index out of bounds.");
    return line_nodes_[i];}

/// For node i {0,1}, give me indices of lines.
template<> template<> inline const IdxVector<1> RefElement<1>::interact<1,0>(unsigned int i)
{   ASSERT(i < RefElement<1>::n_nodes, "Index out of bounds.");
    return node_lines_[i];}

/// For node i {0,1,2}, give me indices of lines.
template<> template<> inline const IdxVector<2> RefElement<2>::interact<1,0>(unsigned int i)
{   ASSERT(i < RefElement<2>::n_nodes, "Index out of bounds.");
    return node_lines_[i];}

/// For node i {0,1,2,3}, give me indices of lines.
template<> template<> inline const IdxVector<3> RefElement<3>::interact<1,0>(unsigned int i)
{   ASSERT(i < RefElement<3>::n_nodes, "Index out of bounds.");
    return node_lines_[i];}
    
/// For side i {0,1,2}, give me indices of its nodes.
template<> template<> inline const IdxVector<3> RefElement<3>::interact<0,2>(unsigned int i)
{   ASSERT(i < RefElement<3>::n_sides, "Index out of bounds.");
    return side_nodes_[i];}

/// For node i {0,1,2,3}, give me indices of sides.
template<> template<> inline const IdxVector<3> RefElement<3>::interact<2,0>(unsigned int i)
{   ASSERT(i < RefElement<3>::n_sides, "Index out of bounds.");
    return node_sides_[i];}
    
/// For line i {0,1,2,3}, give me indices of sides.
template<> template<> inline const IdxVector<2> RefElement<3>::interact<2,1>(unsigned int i)
{   ASSERT(i < RefElement<3>::n_lines, "Index out of bounds.");
    return line_sides_[i];}

/// For side i {0,1,2}, give me indices of its lines.
template<> template<> inline const IdxVector<3> RefElement<3>::interact<1,2>(unsigned int i)
{   ASSERT(i < RefElement<3>::n_sides, "Index out of bounds.");
    return side_lines_[i];}
    
template<unsigned int dim> template<unsigned int OutDim, unsigned int InDim> 
inline const IdxVector< (InDim>OutDim ? InDim+1 : dim-InDim) > RefElement<dim>::interact(unsigned int i)
{   ASSERT_LESS(OutDim, dim);
    ASSERT_LESS(InDim, dim);}

#endif /* REF_ELEMENT_HH_ */
