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
 */

/** Auxilliary class representing vector of indices (unsigned int).
 * @tparam Size is the fixed size of the vector.
 */
/*
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


    static const std::vector< std::vector< std::vector<unsigned int> > > nodes_of_subelements;

	/**
	 * Number of permutations of nodes on sides.
	 * dim   value
	 * -----------
	 * 1     1
	 * 2     2
	 * 3     6
	 */
	static constexpr unsigned int n_side_permutations = (dim+1)*(2*dim*dim-5*dim+6)/6;

	/**
	 * Permutations of nodes on sides.
     * [n_side_permutations][n_nodes_per_side]
	 */
	static const std::vector< std::vector<unsigned int> > side_permutations;

	/**
	 * For a given permutation @p p of nodes finds its index within @p side_permutations.
	 * @param p Permutation of nodes.
	 */
	static unsigned int permutation_index(unsigned int p[n_nodes_per_side]);

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
    static const IdxVector< (InDim>OutDim ? InDim+1 : dim-InDim) > interact( TInteraction<OutDim,InDim> interaction );


private:
    /// Internal part of the interact function.
    template<unsigned int OutDim, unsigned int InDim> 
    static const IdxVector< (InDim>OutDim ? InDim+1 : dim-InDim) > interact_(unsigned int index);
    
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


template<> const IdxVector<1> RefElement<0>::topology_zeros_[];
template<> const IdxVector<2> RefElement<1>::topology_zeros_[];
template<> const IdxVector<3> RefElement<2>::topology_zeros_[];
template<> const IdxVector<6> RefElement<3>::topology_zeros_[];






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
/*
template <unsigned int Size>
IdxVector<Size>::IdxVector(std::array<unsigned int,Size> data_in)
: data_(data_in){}

template <unsigned int Size>
IdxVector<Size>::IdxVector(std::initializer_list<unsigned int> data_in)
{
    ASSERT_EQ_DBG(data_in.size(), Size).error("Incorrect data size.");
    std::copy(data_in.begin(), data_in.end(), data_);
}

template <unsigned int Size>
inline unsigned int IdxVector<Size>::operator[](unsigned int idx) const
{   ASSERT_LT_DBG(idx, Size).error("Index out of bounds.");
    return data_[idx]; }
    
*/

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
template<> template<> inline const IdxVector<1> RefElement<1>::interact_<0,0>(unsigned int i)
{   ASSERT_LT_DBG(i, RefElement<1>::n_nodes).error("Index out of bounds.");
    return IdxVector<1>({i});}

/// For line i {0}, give me indices of its nodes.
template<> template<> inline const IdxVector<2> RefElement<1>::interact_<0,1>(unsigned int i)
{   ASSERT_LT_DBG(i, RefElement<1>::n_lines).error("Index out of bounds.");
    return line_nodes_[i];}

/// For line i {0,1,2}, give me indices of its nodes.
template<> template<> inline const IdxVector<2> RefElement<2>::interact_<0,1>(unsigned int i)
{   ASSERT_LT_DBG(i, RefElement<2>::n_lines).error("Index out of bounds.");
    return line_nodes_[i];}

/// For line i {0,1,2,3,4,5}, give me indices of its nodes.
template<> template<> inline const IdxVector<2> RefElement<3>::interact_<0,1>(unsigned int i)
{   ASSERT_LT_DBG(i, RefElement<3>::n_lines).error("Index out of bounds.");
    return line_nodes_[i];}

/// For node i {0,1}, give me indices of lines.
template<> template<> inline const IdxVector<1> RefElement<1>::interact_<1,0>(unsigned int i)
{   ASSERT_LT_DBG(i, RefElement<1>::n_nodes).error("Index out of bounds.");
    return node_lines_[i];}

/// For node i {0,1,2}, give me indices of lines.
template<> template<> inline const IdxVector<2> RefElement<2>::interact_<1,0>(unsigned int i)
{   ASSERT_LT_DBG(i, RefElement<2>::n_nodes).error("Index out of bounds.");
    return node_lines_[i];}

/// For node i {0,1,2,3}, give me indices of lines.
template<> template<> inline const IdxVector<3> RefElement<3>::interact_<1,0>(unsigned int i)
{   ASSERT_LT_DBG(i, RefElement<3>::n_nodes).error("Index out of bounds.");
    return node_lines_[i];}
    
/// For side i {0,1,2}, give me indices of its nodes.
template<> template<> inline const IdxVector<3> RefElement<3>::interact_<0,2>(unsigned int i)
{   ASSERT_LT_DBG(i, RefElement<3>::n_sides).error("Index out of bounds.");
    return side_nodes_[i];}

/// For node i {0,1,2,3}, give me indices of sides.
template<> template<> inline const IdxVector<3> RefElement<3>::interact_<2,0>(unsigned int i)
{   ASSERT_LT_DBG(i, RefElement<3>::n_sides).error("Index out of bounds.");
    return node_sides_[i];}
    
/// For line i {0,1,2,3}, give me indices of sides.
template<> template<> inline const IdxVector<2> RefElement<3>::interact_<2,1>(unsigned int i)
{   ASSERT_LT_DBG(i, RefElement<3>::n_lines).error("Index out of bounds.");
    return line_sides_[i];}

/// For side i {0,1,2}, give me indices of its lines.
template<> template<> inline const IdxVector<3> RefElement<3>::interact_<1,2>(unsigned int i)
{   ASSERT_LT_DBG(i, RefElement<3>::n_sides).error("Index out of bounds.");
    return side_lines_[i];}
    
template<unsigned int dim> template<unsigned int OutDim, unsigned int InDim> 
inline const IdxVector< (InDim>OutDim ? InDim+1 : dim-InDim) > RefElement<dim>::interact_(unsigned int i)
{
    ASSERT(false)(dim)(OutDim)(InDim)(i).error("Not implemented.");
    //ASSERT_LT_DBG(OutDim, dim);
    //ASSERT_LT_DBG(InDim, dim);
    return IdxVector< (InDim>OutDim ? InDim+1 : dim-InDim) >();  // just to avoid warning for missing return
}


template<unsigned int dim>
template < template <unsigned int OutDim, unsigned int InDim> class TInteraction, unsigned int OutDim, unsigned int InDim>
inline  const IdxVector< (InDim>OutDim ? InDim+1 : dim-InDim) > RefElement<dim>::interact( TInteraction<OutDim,InDim> interaction )
{
    return interact_<OutDim,InDim>(interaction.i_);
}
    
#endif /* REF_ELEMENT_HH_ */
