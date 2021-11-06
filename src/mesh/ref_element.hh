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
 * @file    ref_element.hh
 * @brief   Class RefElement defines numbering of vertices, sides, calculation of normal vectors etc.
 * @author  Jan Stebel
 * @todo
 *
 * TODO:  reconsider following whether it is actual...
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

#include <vector>                      // for vector
#include <array>
#include <armadillo>
#include "system/armor.hh"
#include "system/asserts.hh"


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
 * 0        0                  0        0,1                0        0,1,2      OUT
 * 1        1                  1        0,2                1        0,1,3      IN
 *                             2        1,2                2        0,2,3      OUT
 *                                                         3        1,2,3      IN
 *
 *
 * nodes coordinates:                                                               
 * 0        [0]                0        [0,0]              0        [0,0,0]    
 * 1        [1]                1        [1,0]              1        [1,0,0]    
 *                             2        [0,1]              2        [0,1,0]       
 *                                                         3        [0,0,1]    
 * 
 * barycentric coordinates of nodes:
 * 0        [1,0]              0        [1,0,0]            0        [1,0,0,0]
 * 1        [0,1]              1        [0,1,0]            1        [0,1,0,0]
 *                             2        [0,0,1]            2        [0,0,1,0]
 *                                                         3        [0,0,0,1]
 *
 *
 *
 *  Element node permutation for matching sides
 *
 *  1. 2D case: Can add mesh elements in "layers" so that at most toe edges are prescribed and one is free.
 *     possible cases, two edges A and B considered counter clockwise around a node:
 *     A out, B out : vertex 0 to node, regular orientation (normal up)
 *     A out, B in  : vertex 1 to node, regular orientation
 *     A in , B out : vertex 1 to node, inverted orientation (normal down)
 *     A in , B in  : vertex 2 to node regular orientation
 *
 *     Result: can be matched using two element orientations, element inversion in 1/3 of cases
 *
 *  2. 3d elements around an oriented edge
 *     a. possible face configurations
 *        U: side 0 attached to the edge, positive; face vertices: 0 (edge 0), 1 (edge 1), 2
 *        V: side 1 attached to the edge, negative; face vertices: 0 (edge 0), 2 (edge 1), 1
 *        W: side 2 attached to the edge, positive; face vertices: 1 (edge 0), 2 (edge 1), 0
 *
 *     b. faces of the element edges, configuration of faces attached to the element edge; edge pointing up, looking from outside position, face A on right, face B on left
 *        E0: A (side 0) in  U, B (side 1) in U
 *        E1: A (side 2) in  U, B (side 0) in V
 *        E2: A (side 1) in  V, B (side 2) in V
 *        E3. A (side 0) in  W, B (side 3) in U
 *        E4: A (side 3) in  V, B (side 1) in W
 *        E5: A (side 2) in  W, B (side 3) in W
 *        not presented combination:
 *        A in U, B in W   =  E3 inverted
 *        A in V, B in U   =  E1 inverted
 *        A in W, b in V   =  E4 inverted
 *
 *        Result can permute element veritices to arrange elements around an edge, possibly element inversion in 1/3 of cases
 *
 *  3. 3d elements around a node, up to 3 given faces, one free face
 *     denote V0, V1, V2 vertices of the free face, V3 common vertex of given faces
 *     Denote edges of this element ABCDEF with +- orientation
 *     a. orientation of edges of the vertices
 *     V0: E0+ >V1, E1+ >V2, E2+ >V3    +++
 *     V1: E0- >V0, E3+ >V2, E4+ >V3    -++
 *     V2: E1- >V0, E3- >V1, E5+ >V3    --+
 *     V3: E2- >V0, E4- >V1, E5- >V2    ---
 *
*                   vertices   edges       normal (out = +)
 *   3D - sides: S0: 0 1 2      E0 E1 E3    -
 *               S1: 0 1 3      E0 E2 E4    +
 *               S2: 0 2 3      E1 E2 E5    -
 *               S3: 1 2 3      E3 E4 E5    -
 *
 *        edges: A E0: V0 -> V1  x direction
 *               B E1: V0 -> V2  y direction
 *               C E2: V0 -> V3  z direction
 *               D E3: V1 -> V2
 *               E E4: V1 -> V3
 *               F E5: V2 -> V3

                 CE -> A
                 EF -> D
                 FC -> B
 *
 *   with edges  DEF, other three edges ABC
 *     configuration given by orientation of edges ABC and possibly by orientation of edges DEF
 *     ABC orientation + down, - up; DEF orientation positive if match side 0, i.e. D (AB), E (AC), F (BC)
 *     faces denoted fD,fE,fF, their configurations considered with respect to the edges DEF respectively
 *     front face denoted: ff
 *     vertices V0, V1, V2 top vertices with edges D (V0V1) E (V0V2) F (V1V2); V3 bottom
 *     invalid = corner case, can not continue must either modify one of the faces, or introduce unmatching face
 *     impossible = combination of edges that can not happen
 *     edge cases arranged for the vertex 3
 *
 *     CEF ABD  face config: fD   fE  fF  ff
 *     +++ +++               U    U   U
 *     +++ ++-               U    U   U-
 *     +++ +-+ invalid
 *     +++ +--               U    U   U
 *     +++ -++               U    U   U
 *     +++ -+- invalid              U    U   U-
 *     +++ --+
 *     +++ ---               U    U   U
 *
 *     ++- +++ impossible              U    U   U
 *     ++- ++- impossible              U    U   U-
 *     ++- +-+ impossible
 *     ++- +--               U    U   U
 *     ++- -++ impossible              U    U   U
 *     ++- -+- impossible             U    U   U-
 *     ++- --+ impossible
 *     ++- ---               U    U   U

 *     +-+ +++ impossible              U    U   U
 *     +-+ ++- impossible              U    U   U-
 *     +-+ +-+ impossible
 *     +-+ +--               U    U   U
 *     +-+ -++ impossible              U    U   U
 *     +-+ -+- impossible             U    U   U-
 *     +-+ --+ impossible
 *     +-+ ---               U    U   U

 *     +-- +++ impossible              U    U   U
 *     +-- ++- impossible              U    U   U-
 *     +-- +-+ impossible
 *     +-- +--               U    U   U
 *     +-- -++ impossible              U    U   U
 *     +-- -+- impossible             U    U   U-
 *     +-- --+ impossible
 *     +-- ---               U    U   U

 *     -++ +++ impossible              U    U   U
 *     -++ ++- impossible              U    U   U-
 *     -++ +-+ impossible
 *     -++ +--               U    U   U
 *     -++ -++ impossible              U    U   U
 *     -++ -+- impossible             U    U   U-
 *     -++ --+ impossible
 *     -++ ---               U    U   U

 *     -+- +++ impossible              U    U   U
 *     -+- ++- impossible              U    U   U-
 *     -+- +-+ impossible
 *     -+- +--               U    U   U
 *     -+- -++ impossible              U    U   U
 *     -+- -+- impossible             U    U   U-
 *     -+- --+ impossible
 *     -+- ---               U    U   U

 *     --+ +++ impossible              U    U   U
 *     --+ ++- impossible              U    U   U-
 *     --+ +-+ impossible
 *     --+ +--               U    U   U
 *     --+ -++ impossible              U    U   U
 *     --+ -+- impossible             U    U   U-
 *     --+ --+ impossible
 *     --+ ---               U    U   U

 *     --- +++               U    U   U
 *     --- ++-               U    U   U-
 *     --- +-+ invalid
 *     --- +--               U    U   U
 *     --- -++               U    U   U
 *     --- -+- invalid             U    U   U-
 *     --- --+
 *     --- ---               U    U   U
 *
 *
 *
 */


template<std::size_t Size>
using IdxVector = std::array<unsigned int, Size>;


/** Auxilliary structure that is used to pass template arguments into interact function of RefElement:
 * RefElement<dim>::interact( Interaction<OutDim,InDim>(i) )
 * 
 * This enables automatic deduction of dimensional template arguments.
 * @see @p RefElement<dim>::interact
 */
template <unsigned int OutDim, unsigned int InDim>
struct Interaction {
    Interaction(unsigned int i) : i_(i) {}
    unsigned int i_;
};



class _AuxInteract {
public:
    // Order clockwise, faces opposite to the lines from node_lines.
    // !! dependes on S3 inversion
    static constexpr IdxVector<3> S3_node_sides [2][4]
            = { { { 2, 1, 0 },
                  { 3, 0, 1 },
                  { 3, 2, 0 },
                  { 3, 1, 2 }},
                { { 2, 0, 1 },
                  { 3, 1, 0 },
                  { 3, 0, 2 },
                  { 3, 2, 1 }}};

    // faces adjecent to given edge, first is the right face when looking form outside with the edge pointing up.
    // !! dependes on S3 inversion
    static constexpr IdxVector<2> S3_line_sides [2][6]
        = { { {0,1},
              {2,0},
              {1,2},
              {0,3},
              {3,1},
              {2,3}},
            { {1,0},
              {0,2},
              {2,1},
              {3,0},
              {1,3},
              {3,2}}};

    // Order clockwise looking over the vertex to center; smallest index first
    // !! dependes on S3 inversion
    static constexpr IdxVector<3> S3_node_lines [2][4]
        = { { {0,1,2},
              {0,4,3},
              {1,3,5},
              {2,5,4}},
            { {0,2,1},
              {0,3,4},
              {1,5,3},
              {2,4,5}}};

};


template<unsigned int dim>
class RefElement
{
public:
    typedef arma::vec::fixed<dim> LocalPoint;
    /**
     * Barycentric coordinates.
     *
     * e.g. coordinates (a,b,c) on triangle with vertices X, Y, Z
     * represents a point: a*X+b*Y+c*Z
     */
    typedef Armor::ArmaVec<double, dim+1> BaryPoint;
    typedef Armor::ArmaVec<double, dim> FaceBaryPoint;

    DECLARE_EXCEPTION( ExcInvalidPermutation, << "Side permutation not found.\n" );
        
	/**
	 * Return coordinates of given node.
     * @see the class documentation @p RefElement
	 * @param nid Node number.
     * NOTE: Implementation is dependent on current node and side numbering.
	 */
	static LocalPoint node_coords(unsigned int nid);
    
	/**
	 * Compute normal vector to a given side.
	 * @param sid Side number.
	 */
	static LocalPoint normal_vector(unsigned int sid);


	/**
	 * If the given barycentric coordinate is in the ref. element, return unchanged.
	 * If the given barycentric coordinate is out of the ref. element,
	 * project it on the surface of the ref. element.
	 */
	static BaryPoint clip(const BaryPoint &barycentric);

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
    static const unsigned int n_sides_per_node = dim;       ///< Number of sides with one common line.
    
	/// Number of lines on boundary of one side.
	static const unsigned int n_lines_per_side = (unsigned int)((dim * (dim - 1)) / 2);//( dim == 3 ? 3 : 0);// Kombinační číslo dim nad dvěma

	/// Number of lines, i.e. @p object of dimension @p dim-2 on the boundary of the reference element.
	static const unsigned int n_lines = (unsigned int)((dim * (dim + 1)) / 2); //( dim == 3 ? 6 : dim == 2 ? 3 : dim == 1 ? 1 : 0); součet posloupnosti


    static const std::vector< std::vector< std::vector<unsigned int> > > nodes_of_subelements;


    /** @brief Converts from local to barycentric coordinates.
     * @param lp point in local coordinates (x,y)
     * @return point in barycentric coordinates (1-x-y, x, y)
     */
    static BaryPoint local_to_bary(const LocalPoint& lp);
    
    /** @brief Converts from barycentric to local coordinates.
     * @param bp point in barycentric coordinates
     * @return point in local coordinates
     */
    static LocalPoint bary_to_local(const BaryPoint& bp);
    
	typedef std::vector<BaryPoint> BarycentricUnitVec;

	/**
	 * Used in the clip method.
	 */
	static BarycentricUnitVec make_bary_unit_vec();

    /**
     * For given barycentric coordinates on the ref element returns barycentric coordinates
     * on the ref. element of given face. Assumes that the input point is on the face.
     * Barycentric order: (complanatory, local_coords )
     */
    static FaceBaryPoint barycentric_on_face(const BaryPoint &barycentric, unsigned int i_face);


    typedef const std::vector<LocalPoint> & CentersList;
    static CentersList centers_of_subelements(unsigned int sub_dim);
    
    /**
     * Return (1) number of zeros and (2) positions of zeros in barycentric coordinates.
     * @p tolerance serves for testing zero values of @p barycentric coordinates.
     */
    static std::pair<unsigned int, unsigned int> zeros_positions(const BaryPoint &barycentric,
                                                                 double tolerance = std::numeric_limits<double>::epsilon()*2);
    
    /**
     * According to positions of zeros in barycentric coordinates, it gives the index of subdim-simplex
     * in the reference element. Number of zeros must be equal to (3-subdim).
     * e.g.:
     * if 1 zeros, return index of side (subdim 2)
     * if 2 zeros, return index of edge (subdim 1)
     * if 3 zeros, return index of vertex (subdim 0)
     */
    template<unsigned int subdim> static unsigned int topology_idx(unsigned int zeros_positions);
    
    /** Function returns number of subdim-simplices inside dim-simplex.
     * The aim is covering all the n_**** members with a single function.
     * TODO: think of generalization for n_****_per_**** members, like function @p interact:
     * template<unsigned int subdimA, unsigned int subdimB> static unsigned int count();
     */
    template<unsigned int subdim> static unsigned int count();
    
    /**
     * @param sid - index of a sub-simplex in a simplex
     * return an array of barycentric coordinates on <dim> simplex from <subdim> simplex
     * for example: simplex<3> - ABCD and its subsubsimplex<1> AD (line index: 3)
     * AD has barycoords for A (1,0), for D (0,1), but A in ABCD is (1,0,0,0) and D is (0,0,0,1)
     * this method creates array ((1,0,0,0),(0,0,0,1))
     */
    template<unsigned int subdim> static arma::mat::fixed<dim+1,subdim+1> bary_coords(unsigned int sid);

    /** Interpolate barycentric coords to a higher dimension of a simplex.
     * @param coord - barycentric coords of a point on a sub-simplex
     * @param sub_simplex_idx - id of sub-simplex on a simplex
     */
    template<unsigned int subdim> static BaryPoint interpolate(arma::vec::fixed<subdim+1> coord, int sub_simplex_idx);


    /**
     * Basic line interpolation.
     */
    static BaryPoint line_barycentric_interpolation(BaryPoint first_coords, 
                                                    BaryPoint second_coords, 
                                                    double first_theta, double second_theta, double theta);
    
    /**
     * Usage: 
     * RefElement<3>::interact(Interaction<2,0>(1))
     * (means: In tetrahedron <3>, give indices of sides <2>, connected by node <0> with index 1)
     * RefElement<3>::interact(Interaction<2,0>(1))[1]
     * (as above, but give only the side with index 1)
     * 
     * Template usage: RefElement<dim>::interact(Interaction<OutDim, InDim>(i))[j]
     * (means: on dim-dimensional reference element, go on InDim-dimensional subelement with index i,
     * which connects OutDim-dimnesional subelements and select the one with index j)
     * 
     * This method serves as an interface to topology information of the reference element.
     * It returns indices of OutDim-dimensional object
     * of InDim-dimnesional object of given index
     * in dim-dimnesional reference element.
     * @tparam interaction - auxilliary object carying the index and the template arguments OutDim and InDim
     * @tparam OutDim - output dimension (give me node-0, line-1, side-2), <= dim
     * @tparam InDim - input dimension (for node-0, line-1, side-2), <= dim
     * @return vector of indices of OutDim-dimensional subelements represented by @p IdxVector object.
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
    template < template <unsigned int OutDim, unsigned int InDim> class TInteraction, unsigned int OutDim, unsigned int InDim>
    static const IdxVector< (InDim>OutDim ? InDim+1 : dim-InDim) > interact( TInteraction<OutDim,InDim> interaction, bool inv = false );


private:
    /// Internal part of the interact function.
    template<unsigned int OutDim, unsigned int InDim> 
    static const IdxVector< (InDim>OutDim ? InDim+1 : dim-InDim) > interact_(unsigned int index, bool inv = false);
    
    static const std::vector<IdxVector<n_nodes_per_line>> line_nodes_; ///< [n_lines] For given line, returns its nodes indices.
    static const std::vector<IdxVector<n_lines_per_node>> node_lines_; ///< [n_nodes] For given node, returns lines indices.
    static const std::vector<IdxVector<n_nodes_per_side>> side_nodes_; ///< [n_sides] For given side, returns nodes indices. For @p dim == 3.
    static const std::vector<IdxVector<n_sides_per_node>> node_sides_; ///< [n_nodes] For given node, returns sides indices. For @p dim == 3.
    static const std::vector<IdxVector<n_sides_per_line>> line_sides_; ///< [n_lines] For given line, returns sides indices. For @p dim == 3.
    static const std::vector<IdxVector<n_lines_per_side>> side_lines_; ///< [n_sides] For given side, returns lines indices. For @p dim == 3.

    //TODO: implement for 1d and 2d
    /**
     * Consider an n-face (node, edge, face, bulk) with dimension `subdim` and
     * index within subdimension `idx`. Barycentric coordinates of all points
     * on the n-face have unique pattern of zero coordinates.
     *
     * topology_zeros_[subdim][idx] is a bitfield with '1' where the pattern have zeros.
     */
    static const IdxVector<(n_lines > n_nodes) ? n_lines : n_nodes> topology_zeros_[dim+1];
};


// Declarations of explicit specialization of static memebers.

template<> const IdxVector<1> RefElement<0>::topology_zeros_[];
template<> const IdxVector<2> RefElement<1>::topology_zeros_[];
template<> const IdxVector<3> RefElement<2>::topology_zeros_[];
template<> const IdxVector<6> RefElement<3>::topology_zeros_[];


template<> const std::vector<IdxVector<2>> RefElement<1>::line_nodes_;
template<> const std::vector<IdxVector<2>> RefElement<2>::line_nodes_;
template<> const std::vector<IdxVector<2>> RefElement<3>::line_nodes_;
template<> const std::vector<IdxVector<1>> RefElement<1>::node_lines_;
template<> const std::vector<IdxVector<2>> RefElement<2>::node_lines_;

template<> const std::vector<IdxVector<3>> RefElement<3>::node_lines_;
template<> const std::vector<IdxVector<3>> RefElement<3>::side_nodes_;
template<> const std::vector<IdxVector<3>> RefElement<3>::node_sides_;
template<> const std::vector<IdxVector<2>> RefElement<3>::line_sides_;
template<> const std::vector<IdxVector<3>> RefElement<3>::side_lines_;


template<> const IdxVector<1> RefElement<0>::topology_zeros_[];
template<> const IdxVector<2> RefElement<1>::topology_zeros_[];
template<> const IdxVector<3> RefElement<2>::topology_zeros_[];
template<> const IdxVector<6> RefElement<3>::topology_zeros_[];






// 0: nodes of nodes
// 1: nodes of lines
// 2: nodes of sides
// 3: nodes of tetrahedron
template<> const std::vector< std::vector< std::vector<unsigned int> > > RefElement<0>::nodes_of_subelements;
template<> const std::vector< std::vector< std::vector<unsigned int> > > RefElement<1>::nodes_of_subelements;
template<> const std::vector< std::vector< std::vector<unsigned int> > > RefElement<2>::nodes_of_subelements;
template<> const std::vector< std::vector< std::vector<unsigned int> > > RefElement<3>::nodes_of_subelements;









/************************* template implementation ****************************/

template<unsigned int dim>
template<unsigned int subdim> inline
arma::mat::fixed<dim+1,subdim+1> RefElement<dim>::bary_coords(unsigned int sid){
        ASSERT_LT_DBG(subdim, dim).error("Dimension mismatch!");
        arma::mat::fixed<dim+1,subdim+1> bary_c;
        
        for(unsigned int i = 0; i < subdim+1; i++){
        	unsigned int nid = interact_<0,subdim>(sid)[i];
            bary_c.col(i).zeros();
            bary_c.col(i)(nid) = 1;
        }       
    
        return bary_c;
}


template<unsigned int dim> inline
arma::vec::fixed<dim> RefElement<dim>::node_coords(unsigned int nid)
{
	ASSERT_LT_DBG(nid, n_nodes).error("Node number is out of range!");

	arma::vec::fixed<dim> p;
	p.zeros();

    if (nid > 0)
        p(nid-1) = 1;

	return p;
}


template<unsigned int dim>
template<unsigned int subdim> 
auto RefElement<dim>::interpolate(arma::vec::fixed<subdim+1> coord, int sub_simplex_idx) -> BaryPoint
{
    return RefElement<dim>::bary_coords<subdim>(sub_simplex_idx)*coord;
}


template<> template<> inline unsigned int RefElement<3>::count<0>()
{ return n_nodes; }
template<> template<> inline unsigned int RefElement<3>::count<1>()
{ return n_lines; }
template<> template<> inline unsigned int RefElement<3>::count<2>()
{ return n_sides; }
template<> template<> inline unsigned int RefElement<3>::count<3>()
{ return 1; }
template<> template<> inline unsigned int RefElement<2>::count<0>()
{ return n_nodes; }
template<> template<> inline unsigned int RefElement<2>::count<1>()
{ return n_lines; }
template<> template<> inline unsigned int RefElement<2>::count<2>()
{ return 1; }
template<> template<> inline unsigned int RefElement<2>::count<3>()
{ return 0; }
template<> template<> inline unsigned int RefElement<1>::count<0>()
{ return n_nodes; }
template<> template<> inline unsigned int RefElement<1>::count<1>()
{ return 1; }
template<> template<> inline unsigned int RefElement<1>::count<2>()
{ return 0; }
template<> template<> inline unsigned int RefElement<1>::count<3>()
{ return 0; }
template<> template<> inline unsigned int RefElement<0>::count<0>()
{ return 1; }
template<> template<> inline unsigned int RefElement<0>::count<1>()
{ return 0; }
template<> template<> inline unsigned int RefElement<0>::count<2>()
{ return 0; }
template<> template<> inline unsigned int RefElement<0>::count<3>()
{ return 0; }

template<unsigned int dim>
template<unsigned int subdim>
unsigned int RefElement<dim>::topology_idx(unsigned int zeros_positions)
{
    for(unsigned int i=0; i < RefElement<dim>::count<subdim>(); i++){
        if(zeros_positions == topology_zeros_[subdim][i]) return i;
    }
    ASSERT(0).error("Undefined zero pattern.");
    return -1;
}


/// This function is for "side_nodes" - for given side, give me nodes (0->0, 1->1).
template<> template<> inline const IdxVector<1> RefElement<1>::interact_<0,0>(unsigned int i, bool inv)
{   ASSERT_LT_DBG(i, RefElement<1>::n_nodes)(inv).error("Index out of bounds.");
    return IdxVector<1>({i});}

/// For line i {0}, give me indices of its nodes.
template<> template<> inline const IdxVector<2> RefElement<1>::interact_<0,1>(unsigned int i, bool inv)
{   ASSERT_LT_DBG(i, RefElement<1>::n_lines)(inv).error("Index out of bounds.");
    return line_nodes_[i];}

/// For line i {0,1,2}, give me indices of its nodes.
template<> template<> inline const IdxVector<2> RefElement<2>::interact_<0,1>(unsigned int i, bool inv)
{   ASSERT_LT_DBG(i, RefElement<2>::n_lines)(inv).error("Index out of bounds.");
    return line_nodes_[i];}

/// For line i {0,1,2,3,4,5}, give me indices of its nodes.
template<> template<> inline const IdxVector<2> RefElement<3>::interact_<0,1>(unsigned int i, bool inv)
{   ASSERT_LT_DBG(i, RefElement<3>::n_lines)(inv).error("Index out of bounds.");
    return line_nodes_[i];}

/// For node i {0,1}, give me indices of lines.
template<> template<> inline const IdxVector<1> RefElement<1>::interact_<1,0>(unsigned int i, bool inv)
{   ASSERT_LT_DBG(i, RefElement<1>::n_nodes)(inv).error("Index out of bounds.");
    return node_lines_[i];}

/// For node i {0,1,2}, give me indices of lines.
template<> template<> inline const IdxVector<2> RefElement<2>::interact_<1,0>(unsigned int i, bool inv)
{   ASSERT_LT_DBG(i, RefElement<2>::n_nodes)(inv).error("Index out of bounds.");
    return node_lines_[i];}

/// For node i {0,1,2,3}, give me indices of lines.
template<> template<> inline const IdxVector<3> RefElement<3>::interact_<1,0>(unsigned int i, bool inv)
{   ASSERT_LT_DBG(i, RefElement<3>::n_nodes).error("Index out of bounds.");
    return _AuxInteract::S3_node_lines[inv][i];}
    
/// For side i {0,1,2}, give me indices of its nodes.
template<> template<> inline const IdxVector<3> RefElement<3>::interact_<0,2>(unsigned int i, bool inv)
{   ASSERT_LT_DBG(i, RefElement<3>::n_sides)(inv).error("Index out of bounds.");
    return side_nodes_[i];}

/// For node i {0,1,2,3}, give me indices of sides.
template<> template<> inline const IdxVector<3> RefElement<3>::interact_<2,0>(unsigned int i, bool inv)
{   ASSERT_LT_DBG(i, RefElement<3>::n_sides).error("Index out of bounds.");
    return _AuxInteract::S3_node_sides[inv][i];}
    
/// For line i {0,1,2,3}, give me indices of sides.
template<> template<> inline const IdxVector<2> RefElement<3>::interact_<2,1>(unsigned int i, bool inv)
{   ASSERT_LT_DBG(i, RefElement<3>::n_lines).error("Index out of bounds.");
    return _AuxInteract::S3_line_sides[inv][i];}

/// For side i {0,1,2}, give me indices of its lines.
template<> template<> inline const IdxVector<3> RefElement<3>::interact_<1,2>(unsigned int i, bool inv)
{   ASSERT_LT_DBG(i, RefElement<3>::n_sides)(inv).error("Index out of bounds.");
    return side_lines_[i];}
    
template<unsigned int dim> template<unsigned int OutDim, unsigned int InDim> 
inline const IdxVector< (InDim>OutDim ? InDim+1 : dim-InDim) > RefElement<dim>::interact_(unsigned int i, bool inv)
{
    ASSERT(false)(dim)(OutDim)(InDim)(i)(inv).error("Not implemented.");
    //ASSERT_LT_DBG(OutDim, dim);
    //ASSERT_LT_DBG(InDim, dim);
    return IdxVector< (InDim>OutDim ? InDim+1 : dim-InDim) >();  // just to avoid warning for missing return
}


template<unsigned int dim>
template < template <unsigned int OutDim, unsigned int InDim> class TInteraction, unsigned int OutDim, unsigned int InDim>
inline  const IdxVector< (InDim>OutDim ? InDim+1 : dim-InDim) > RefElement<dim>::interact( TInteraction<OutDim,InDim> interaction , bool inv)
{
    return interact_<OutDim,InDim>(interaction.i_, inv);
}
    
#endif /* REF_ELEMENT_HH_ */
