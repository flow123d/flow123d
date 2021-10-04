/*!
 *
﻿* Copyright (C) 2015 Technical University of Liberec.  All rights reserved.
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
 * @file    intersection_point_aux.hh
 * @brief   Internal class representing intersection point.
 * @author  Viktor Fris, Pavel Exner
 *
 */

#ifndef INTERSECTIONPOINT_H_
#define INTERSECTIONPOINT_H_

#include <armadillo>
template <int spacedim> class ElementAccessor;



static const double geometry_epsilon = 1e-9;

//http://en.cppreference.com/w/cpp/types/numeric_limits/epsilon
// static const double rounding_epsilon = 8*std::numeric_limits<double>::epsilon();

//forward declare
template<unsigned int N, unsigned int M> class IntersectionPointAux;
template<unsigned int N, unsigned int M> std::ostream& operator<<(std::ostream& os, const IntersectionPointAux<N,M>& IP);
    

enum class IntersectionResult {
    positive = 1,
    negative = 0,
    degenerate = 2,
    none = 3
};

inline bool operator<(IntersectionResult a, IntersectionResult b) {
    return int(a) < int(b);
}

/** @brief Internal auxiliary class represents an intersection point of simplex<N> and simplex<M>.
 * 
 * It contains barycentric coordinates of the point on both simplices.
 * Further, it contains topology information of the intersection point.
 * Namely its orientation (according to Plucker coordinates product),
 * pathologic flag, local indices and dimension of element objects (node/line/triangle).
 * Refered as 'IP' in further documentation.
 */
template<unsigned int N, unsigned int M>
class IntersectionPointAux {
    
	arma::vec::fixed<N+1> local_bcoords_A_; ///< Barycentric coordinates of an IP on simplex<N>.
	arma::vec::fixed<M+1> local_bcoords_B_; ///< Barycentric coordinates of an IP on simplex<M>.

	/** @brief Local indices of element objects that intersects.
     * 
     * Can be indices of point (vertex), lines (edge), triangle (side) in intersecting element,
     * depending on the dimensions \p dim_A_, resp. \p dim_B_.
     * 
     * Examples: 
     * dim_A=0, N=2 -> idx_A is local node index of triangle
     * dim_B=2, N=3 -> idx_B is local side index of tetrahedron
     * dim_B=1, N=3 -> idx_B is local line index of tetrahedron
     */
	//@{
	unsigned int idx_A_;
	unsigned int idx_B_;
    //@}

    /** @brief Intersection result according to Plucker products.
     * 
     * All Plucker products > 0 then @b positive. 
     * All Plucker products < 0 then @b negative.
     * @b degenerate case.
     * @b none - no intersection.
     */
	IntersectionResult result_;
	
	unsigned int dim_A_;    ///< Dimension of the object A of intersection. Equal \p N, by default.
    unsigned int dim_B_;    ///< Dimension of the object B of intersection. Equal \p M, by default.
    
public:

    IntersectionPointAux();    ///< Default constructor.
    ~IntersectionPointAux(){}; ///< Destructor.
    
    /** @brief Constructor taking barycentric coordinates on simplices as input parameters.
     * @param lcA barycentric coordinates of IP in Simplex<N>
     * @param lcB barycentric coordinates of IP in Simplex<M>
     * @param dim_A dimension of object A
     * @param dim_B dimension of object B
     */
	IntersectionPointAux(const arma::vec::fixed<N+1> &lcA, const arma::vec::fixed<M+1> &lcB,
                      unsigned int dim_A = N, unsigned int dim_B = M);

	/** @brief Constructor interpolates the second bary coords 
     * of IntersectionPointAux<N,M-1> to IntersectionPointAux<N,M>
     * 
     * Allowed only from dimension \p M 1 to 2 and from 2 to 3.
     * @param IP intersection point of lower dimension of object B
     * @param idx_B is the index of object B of IntersectionPointAux<N,M-1> in object B of IntersectionPointAux<N,M>
     */
	IntersectionPointAux(const IntersectionPointAux<N,M-1> &IP, unsigned int idx_B);

	/** @brief Constructor interpolates the second bary coords
     * of IntersectionPointAux<N,M-2> to IntersectionPointAux<N,M>
     * 
	 * Allowed only from dimension \p M 1 to 3.
     * @param IP intersection point of lower dimension of object B
     * @param idx_B is the index of object B of IntersectionPointAux<N,M-2> in object B of IntersectionPointAux<N,M>
	 */
	IntersectionPointAux(const IntersectionPointAux<N,M-2> &IP, unsigned int idx_B);

    /// Resets the object to default values.
    void clear();
    
    /// Switches the object A and B.
    IntersectionPointAux<M, N> switch_objects() const;
    
    ///@name Setters.
    //@{
    /** @brief Setter for coordinates.
     * @param lcA braycentric coordinates of A
     * @param lcB braycentric coordinates of B
     */
    void set_coordinates(const arma::vec::fixed<N+1> &lcA, const arma::vec::fixed<M+1> &lcB);
    
    /// Setter for topology data. See description of @p idx_A_ and @p dim_A_.
    void set_topology(unsigned int idx_A, unsigned int dim_A,
                      unsigned int idx_B, unsigned int dim_B);
    
    void set_topology_A(unsigned int idx, unsigned int dim_A);  ///< Sets the topology of object A in Simplex<N>.
    void set_topology_B(unsigned int idx, unsigned int dim_B);  ///< Sets the topology of object B in Simplex<M>.
    void set_result(IntersectionResult result);            ///< Setter orientation flag.
    
    //@}


    ///@name Getters.
    //@{
    /// Returns barycentric coordinates in the Simplex<N>.
    const arma::vec::fixed<N+1> &local_bcoords_A() const;
    
    /// Returns barycentric coordinates in the Simplex<M>.
    const arma::vec::fixed<M+1> &local_bcoords_B() const;

    /**
     * The dimension of the X-face the IP lies on within the A Simplex<N>. E.g. IP on and edge of tetrahedra (N=3)
     * have dim_A = 1.
     */
    unsigned int dim_A() const;     ///< Returns dimension of object A.
    /**
     * Same for the B Simplex<M>.
     */
    unsigned int dim_B() const;     ///< Returns dimension of object B.
    /**
     * Index of the X-face of the simplex A.
     */
    unsigned int idx_A() const;     ///<  Returns the index of Simplex<N>.
    /**
     * Index of the X-face of the simplex B.
     */
    unsigned int idx_B() const;     ///<  Returns the index of Simplex<M>.

    /// Result: 0 - negative sign, 1 - positive sign, 2 - degenerate (zero for all sides), 3 - none
    IntersectionResult result() const;   ///<  Returns the orientation.
    //@}
    
    /// Computes real coordinates of IP, given the element @p ele in which IP lies.
    arma::vec::fixed<3> coords(ElementAccessor<3> ele) const;
    
    /// Returns true, if @p other intersection point has the same topology.
    bool topology_equal(const IntersectionPointAux<N,M> &other) const;
    
	/** @brief Comparison operator for sorting the IPs in convex hull tracing algorithm.
     * Compares the points by x-coordinate (in case of a tie, compares by y-coordinate).
	 */
	//bool operator<(const IntersectionPointAux<N,M> &ip) const;

    /// Friend output operator.
    friend std::ostream& operator<< <>(std::ostream& os, const IntersectionPointAux<N,M>& IP);
};


/********************************************* IMPLEMENTATION ***********************************************/

template<unsigned int N, unsigned int M>
void IntersectionPointAux<N,M>::set_coordinates(const arma::vec::fixed< N + 1  >& lcA, const arma::vec::fixed< M + 1  >& lcB)
{   local_bcoords_A_ = lcA;
    local_bcoords_B_ = lcB; }

template<unsigned int N, unsigned int M>
void IntersectionPointAux<N,M>::set_topology(unsigned int idx_A,unsigned int dim_A, unsigned int idx_B, unsigned int dim_B)
{   idx_A_ = idx_A;
    idx_B_ = idx_B;
    dim_A_ = dim_A;
    dim_B_ = dim_B;
}
    
template<unsigned int N, unsigned int M>
void IntersectionPointAux<N,M>::set_result(IntersectionResult res)
{   result_ = res; }

template<unsigned int N, unsigned int M>
unsigned int IntersectionPointAux<N,M>::dim_A() const
{   return dim_A_; }

template<unsigned int N, unsigned int M>
unsigned int IntersectionPointAux<N,M>::dim_B() const
{   return dim_B_; }

template<unsigned int N, unsigned int M>
const arma::vec::fixed< N + 1  >& IntersectionPointAux<N,M>::local_bcoords_A() const
{   return local_bcoords_A_; }

template<unsigned int N, unsigned int M>
const arma::vec::fixed< M + 1  >& IntersectionPointAux<N,M>::local_bcoords_B() const
{   return local_bcoords_B_; }

template<unsigned int N, unsigned int M>
void IntersectionPointAux<N,M>::set_topology_A(unsigned int idx_A, unsigned int dim_A)
{   idx_A_ = idx_A;
    dim_A_ = dim_A; }

template<unsigned int N, unsigned int M>
void IntersectionPointAux<N,M>::set_topology_B(unsigned int idx_B, unsigned int dim_B)
{   idx_B_ = idx_B;
    dim_B_ = dim_B; }

template<unsigned int N, unsigned int M>
unsigned int IntersectionPointAux<N,M>::idx_A() const
{   return idx_A_; }

template<unsigned int N, unsigned int M>
unsigned int IntersectionPointAux<N,M>::idx_B() const
{   return idx_B_; }

template<unsigned int N, unsigned int M>
IntersectionResult IntersectionPointAux<N,M>::result() const
{   return result_; }

#endif /* INTERSECTIONPOINT_H_ */
