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
 * @file    compute_intersection.hh
 * @brief   Fundamental simplicial intersections.
 * @author  Viktor Fris, Pavel Exner
 *
 * Contains classes that compute the fundamental intersection of simplical objects:
 * 1D-2D, 2D-2D, 1D-3D, 2D-3D
 * 
 *
 * TODO: Idea: constructor creates empty object; for each abscissa and triangle only update data 
 * [abscissa, triangle] and (plucker coordinates, products) if not final dimension.
 * Remove init() methods.
 * Wouldn't a flag for final dimension object help code structure?
 * 
 * TODO: 2D-2D and 2D-3D (creating final objects) seems to not need functions for sequence:
 * ComputeIntersection() -remove default
 * init() - make private, call in constructor
 * set_data() - remove, use only proper constructor
 * compute () - can be called in constructor, but keep it in the same fashion as in other classes
 * 
 */

#ifndef COMPUTE_INTERSECTION_H_
#define COMPUTE_INTERSECTION_H_

#include "system/system.hh"
#include "mesh/ref_element.hh"
#include "intersection/intersection_point_aux.hh"


// forward declare
template <int spacedim> class ElementAccessor;
template<unsigned int, unsigned int> class ComputeIntersection;
class Plucker;
template<unsigned int, unsigned int> class IntersectionAux;
template<unsigned int, unsigned int> class IntersectionPointAux;

/// Auxiliary value for Plucker product. If equal this value, it is supposed not to be set yet.
static const double plucker_empty = std::numeric_limits<double>::infinity();

//////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                      CLASS  FOR  1D - 2D
//////////////////////////////////////////////////////////////////////////////////////////////////////////////

/**
 * 
 * @brief Class for 1D-2D intersections.
 * 
 * Computes the intersection of an abscissa and a triangle.
 * Can be used both as sub-algorithm to higher dimnesional intersection
 * and also for 1D-2D intersection as a result.
 */
template<> class ComputeIntersection<1,2> {
public:
    typedef IntersectionPointAux<1,2> IPAux12;
    
    /// Default constructor. Use when this is NOT final intersection object.
	ComputeIntersection();
    /** @brief Constructor, sets abscissa and triangle object.
     * Use when this is final intersection object.
     * It allocates memory, computes plucker coordinates and products.
     */
    ComputeIntersection(ElementAccessor<3> abscissa, ElementAccessor<3> triangle);
	~ComputeIntersection();
    
    /** @brief Computes intersection points of line and triangle.
     * 
     * If Plucker products are non-zero and with the same sign, then IP is inside the triangle.
     * If some of the Plucker products are zero:
     * 1 zero product -> IP is on the triangle side
     * 2 zero products -> IP is at the vertex of triangle (there is no other IP)
     * 3 zero products: -> line and triangle are coplanar
     * 
     * IP is intersection of triangle and whole line (bisector).
     * @param IP12s - input/output vector of IPs. If IP found, it is pushed back. Should be empty on start.
     * @return Orientation flag (0,1 sign of product if get intersection, 2 - three zero products (degenerated),
     * 3 - no intersection
     * 
     * NOTE: Why this is not done in constructor?
     * Because default constructor is called in 1d-3d, 2d-3d and compute() is called later.
     */
	IntersectionResult compute(IPAux12 &IP);
    
    /** Computes final 1d-2d intersection. (Use when this is the resulting dimension object).
     * @param IP12s input/output vector of IPs. If IP found, it is pushed back. Should be empty on start.
     * @return number of intersection points found
     */
    unsigned int compute_final(std::vector<IPAux12> &IP12s);
    
    /** Computes final 1d-2d intersection, supposing situation in 2d plane (only degenerate case).
     * (Use when this is the resulting dimension object).
     * @param IP12s input/output vector of IPs. If IP found, it is pushed back. Should be empty on start.
     * @return number of intersection points found
     */
    unsigned int compute_final_in_plane(std::vector<IPAux12> &IP12s);
    
    /// @name Setters and Getters
    //@{ 

    /// Sets the pointer to Plucker coordinates of the abscissa.
    void set_pc_abscissa(Plucker *p){
		plucker_coordinates_abscissa_ = p;
	}
	
	/** @brief Sets the pointer to Plucker coordinates of the triangle side.
     * @param p pointer to Plucker coordinates
     * @param side_idx local index of a side of triangle
     */
    void set_pc_triangle(Plucker *p, unsigned int side_idx){
		plucker_coordinates_triangle_[side_idx] = p;
	}

	/// Gets the pointer to Plucker coordinates of the abscissa.
    Plucker *get_pc_abscissa(){
		return plucker_coordinates_abscissa_;
	}

	/** @brief Gets the pointer to Plucker coordinates of the triangle side.
	 * @param side_idx local index of a side of triangle
	 */
    Plucker *get_pc_triangle(unsigned int side_idx){
		return plucker_coordinates_triangle_[side_idx];
	}

	/** @brief Sets the Plucker product of abscissa and triangle side (use if computed before).
     * @param value value of the Plucker product
     * @param side_idx local index of a side of triangle
     */
    void set_plucker_product(double* number, unsigned int side_idx)
    { plucker_products_[side_idx] = number; };
    
    /** @brief Getter for Plucker product of abscissa and triangle side.
     * @param side_idx local index of a side of triangle
     * @return pointer to Plucker product
     */
    double* get_plucker_product(unsigned int side_idx)
    { return plucker_products_[side_idx]; };

    /// Gets true if the intersection has been computed already (e.g. in case of IP in vertex).
    bool is_computed() { return computed_; }
    /// Sets the 'computed' flag true. Means that intersection has been computed already (e.g. in case of IP in vertex).
    void set_computed() { computed_ = true; }
    
	//@}

	/** Checks the position of IP on abscissa, possibly changes topology to end points. 
     * Returns -2, -1, 0, 1, 2; -1.
     * 0 - inside abscissa
     * -2, 2 - outside of abscissa
     * -1, 1 - end points
     */
	int check_abscissa_topology(IPAux12& IP);
	
    /// Prints out all the Plucker coordinates.
	void print_plucker_coordinates(std::ostream &os);

    /// Resets Plucker products to 'nullptr'.
    /// Use this CAREFULLY as it does not destroy the objects.
    /// It is intended to be used only from higher dimensions when before destroying.
    void clear_all();
    
private:
    /** Computes Plucker coordinates (abscissa, triangle lines) and Plucker products, 
     * if not computed already, or set from outside.
     */
    void compute_plucker_products();
    
    /** Auxiliary function providing correctly signed Plucker product
     * according to the triangle side orientation.
     * @param i local index of a side of triangle.
     */
    double signed_plucker_product(unsigned int i)
    { return RefElement<2>::normal_orientation(i) ? -(*plucker_products_[i]) : (*plucker_products_[i]); }
    
    /** @brief Computes intersection when nonezero Plucker products are of the same sign.
     * @param IP intersection point (if found)
     * @param local local coordinates of IP (got from Plucker products)
     * @return true, if intersection is found; false otherwise
     */
    bool compute_plucker(IPAux12 &IP, const arma::vec3 &local);
    
    /** Computes intersection of abscissa and triangle side in degenerate case (all products are zero).
     * Inspects also topology.
     * @param side_idx is the local index of the triangle side
     * @param IP is the intersection point (if found)
     * @return true, if intersection is found; false otherwise
     */
    bool compute_degenerate(unsigned int side_idx, IPAux12 &IP);
    
    /** @brief After interpolation, the topology information in tetrahedron must be updated.
     * @param ip intersection point to be corrected
     */
    void correct_triangle_ip_topology(double t, unsigned int i, std::vector<IPAux12> &ip);
    
    /// Flag 'computed'; is true if intersection has been computed already.
    bool computed_;
    ///
    double scale_line_, scale_triangle_;

    /// Pointer to plucker coordinates of abscissa.
	Plucker* plucker_coordinates_abscissa_;
    /// Vector of pointers to plucker coordinates of triangle sides.
	std::vector<Plucker *> plucker_coordinates_triangle_;
    /// Pointers to Plucker products of abscissa and triangle side.
	std::vector<double *> plucker_products_;

};



//////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                      CLASS  FOR  2D - 2D
//////////////////////////////////////////////////////////////////////////////////////////////////////////////

/**
 * @brief Class for 2D-2D intersections.
 * 
 * Computes the intersection of triangle A and triangle B.
 * Uses 6 ComputeIntersection<1,2>:
 *  - 3x side of triangle A X triangle B
 *  - 3x side of triangle B X triangle A
 * 
 * It creates a final intersection.
 * 
 * TODO: do not remember, if the function correct_triangle_ip_topology used to have some meaning here..
 */
template<> class ComputeIntersection<2,2> {
public:
    typedef IntersectionPointAux<1,2> IPAux12;
    typedef IntersectionPointAux<2,2> IPAux22;
    
    /** @brief Default constructor, creates empty object.
     * Resizes vectors for Plucker coordinates and products.
     */
    ComputeIntersection();
    /** @brief Constructor, sets both triangle objects.
     * It allocates memory, computes plucker coordinates and products.
     */
    ComputeIntersection(ElementAccessor<3> triaA, ElementAccessor<3> triaB);
    ~ComputeIntersection();
    
    /** @brief Initializes lower dimensional objects.
     * Sets correctly the pointers to Plucker coordinates and products.
     */
    void init();
    
    /** @brief Computes final 2D-2D intersection.
     * Computes CIs (side vs triangle) in following order: [A0_B, A1_B, A2_B, B0_A, B1_A, B2_A].
     * During cycle determine topology information and set computed flags accordingly.
     * @param intersection final 2D-2D intersection object (output)
     * @return number of intersection points found
     */
    unsigned int compute(IntersectionAux<2,2> &intersection);
    
     /// @name Setters and Getters
    //@{
    
    /** Sets the pointer to Plucker coordinates of the triangle A side.
     * @param p pointer to Plucker coordinates
     * @param side_idx local index of a side of triangle A
     */
    void set_pc_triaA(Plucker *p, unsigned int side_idx){
        plucker_coordinates_[side_idx] = p;
    }
    
    /** Sets the pointer to Plucker coordinates of the triangle B side.
     * @param p pointer to Plucker coordinates
     * @param side_idx local index of a side of triangle B
     */
    void set_pc_triaB(Plucker *p, unsigned int side_idx){
        plucker_coordinates_[3 + side_idx] = p;
    }

    /// Gets the pointer to Plucker coordinates of the triangle A side of given @p side_idx.
    Plucker *get_pc_triaA(unsigned int side_idx){
        return plucker_coordinates_[side_idx];
    }
    
    /// Gets the pointer to Plucker coordinates of the triangle B side of given @p side_idx.
    Plucker *get_pc_triaB(unsigned int side_idx){
        return plucker_coordinates_[3 + side_idx];
    }
    
    /** Sets the pointer to Plucker product of triangle A side and triangle B side.
     * @param number pointer to value of Plucker product
     * @param sideA_idx local index of a side of triangle A
     * @param sideB_idx local index of a side of triangle B
     */
    void set_plucker_product(double * number, unsigned sideA_idx, unsigned int sideB_idx){
        plucker_products_[sideA_idx*3 + sideB_idx] = number;
    }
    /// Gets the pointer to Plucker product of triangle A side @p sideA_idx and triangle B side @p sideB_idx.
    double* get_plucker_product(unsigned sideA_idx, unsigned int sideB_idx){
        return plucker_products_[sideA_idx*3 + sideB_idx];
    }
    //@}
    
    /// Prints out the Plucker coordinates of abscissa and tetrahedron edges.
    void print_plucker_coordinates(std::ostream &os);
    /// Prints out the Plucker coordinates of tetrahedron edges in a tree of tetrahedron sides (triangles).
    void print_plucker_coordinates_tree(std::ostream &os);

    /// Resets Plucker products to 'nullptr'.
    /// Use this CAREFULLY as it does not destroy the objects.
    /// It is intended to be used only from higher dimensions when before destroying.
    void clear_all();
    
private:
    /// Pointers to plucker coordinates of sides of both triangles [triaA[3], triaB[3]], size 6.
    std::vector<Plucker *> plucker_coordinates_;
    /// Pointers to Plucker products of triangles sides [3x[sideA x triaB]]], size 9.
    std::vector<double *> plucker_products_;
    /// Compute intersection for side x triangle [3x[sideA x tria B],3x[sideB x triaA]].
    ComputeIntersection<1,2> CI12[6];
    
    // After interpolation, the topology information in triangle must be updated.
//     void correct_triangle_ip_topology(IntersectionPointAux<2,2> &ip);
};


//////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                      CLASS  FOR  1D - 3D
//////////////////////////////////////////////////////////////////////////////////////////////////////////////

/**
 * @brief * @brief Class for 1D-3D intersections.
 *
 * Computes the intersection of an abscissa and a tetrahedron.
 * Uses 4 ComputeIntersection<1,2> for abscissa and tetrahedron sides intersections.
 */
template<> class ComputeIntersection<1,3> {

public:
    typedef IntersectionPointAux<1,3> IPAux;

	/// Default constructor. Use when this is NOT final intersection object.
    ComputeIntersection();
    /** @brief Constructor, sets abscissa and tetrahedron object.
     * Use when this is final intersection object.
     * It allocates memory, computes plucker coordinates and products.
     */
    ComputeIntersection(ElementAccessor<3> abscissa, ElementAccessor<3> tetrahedron);
	~ComputeIntersection();
	
    /** @brief Initializes lower dimensional objects.
     * Sets correctly the pointers to Plucker coordinates and products.
     */
	void init();
    
    /** @brief Computes intersection points for 1D-3D intersection.
     * Computes lower dimensional CIs abscissa vs tetrahedron face.
     * During cycle determine topology information and set computed flags accordingly.
     * We can have 2 IPs at maximum.
     * If there are 2 IPs:
     * - sort them in direction of abscissa
     * - possibly cut (move IPs so they are inside abscissa)
     * - after cutting interpolate 3D coordinates
     * - check that IPs are not the same, if true, throw one away
     * @param IP13d vector of intersection points (output)
     * @return number of intersection points found
     */
    unsigned int compute(std::vector<IPAux> &IP13s);
    
    /** @brief Computes final 1D-3D intersection.
     * Computes IPs and check if any of them are pathologic to set the resulting object also pathologic.
     * Calls @p compute function inside.
     * @param intersection final 1D-3D intersection object (output)
     * @return number of intersection points found
     */
    unsigned int compute(IntersectionAux<1,3> &intersection);
    
     /// @name Setters and Getters
    //@{ 
    
    /// Sets the pointer to Plucker coordinates of the abscissa.
    void set_pc_abscissa(Plucker *p){
		plucker_coordinates_abscissa_ = p;
	}
	
	/** @brief Sets the pointer to Plucker coordinates of the tetrahedron edge.
	 * @param p pointer to Plucker coordinates
     * @param edge_idx local index of an edge of tetrahedron
     */
    void set_pc_tetrahedron(Plucker *p, unsigned int edge_idx){
		plucker_coordinates_tetrahedron[edge_idx] = p;
	}

	/// Gets the pointer to Plucker coordinates of the abscissa.
    Plucker *get_pc_abscissa(){
		return plucker_coordinates_abscissa_;
	}
    
    /// Gets the pointer to Plucker coordinates of the tetrahedron edge of given @p edge_idx.
    Plucker *get_pc_tetrahedron(unsigned int edge_idx){
		return plucker_coordinates_tetrahedron[edge_idx];
    }
    
    /** @brief Sets the pointer to Plucker product of abscissa and tetrahedron edge.
     * @param number pointer to value of Plucker product
     * @param edge_idx local index of an edge of tetrahedron
     */
    void set_plucker_product(double * number, unsigned edge_idx){
        plucker_products_[edge_idx] = number;
    }
    /// Gets the pointer to Plucker product of abscissa and tetrahedron edge of given @p edge_idx.
    double* get_plucker_product(unsigned edge_idx){
        return plucker_products_[edge_idx];
    }
    //@}
    
	/// Prints out the Plucker coordinates of abscissa and tetrahedron edges.
	void print_plucker_coordinates(std::ostream &os);
    /// Prints out the Plucker coordinates of tetrahedron edges in a tree of tetrahedron sides (triangles).
	void print_plucker_coordinates_tree(std::ostream &os);

    /// Resets Plucker products to 'nullptr'.
    /// Use this CAREFULLY as it does not destroy the objects.
    /// It is intended to be used only from higher dimensions when before destroying.
    void clear_all();
    
private:
    /** @brief After interpolation, the topology information in tetrahedron must be updated.
     * @param ip intersection point to be corrected
     */
    void correct_tetrahedron_ip_topology(double t, unsigned int i, std::vector<IPAux> &ip);
    
    /// Pointer to plucker coordinates of abscissa.
    Plucker* plucker_coordinates_abscissa_;
    /// Vector of pointers to plucker coordinates of tetrahedron edges.
    std::vector<Plucker *> plucker_coordinates_tetrahedron;
    /// Pointers to Plucker products of abscissa and tetrahedron edges.
    std::vector<double *> plucker_products_;
    /**
     * Intersection objects for the tetrahedron's faces.
     * Faces follows element node ordering, element inversion not handled at this stage.
     * Element inversin is taken into account in the compute method.
     * TODO: reuse these calculactions for elements with common face.
     */
    ComputeIntersection<1,2> CI12[4];
};

//////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                      CLASS  FOR  2D - 3D
//////////////////////////////////////////////////////////////////////////////////////////////////////////////

/**
 * @brief Class for 2D-2D intersections.
 * 
 * Computes the intersection of a triangle and a tetrahedron.
 * Uses 3 ComputeIntersection<1,3> for triangle sides vs tetrahedron intersections.
 * Uses 6 ComputeIntersection<1,2> for tetrahedron sides vs triangle intersections.
 */
template<> class ComputeIntersection<2,3> {
public:
    typedef IntersectionPointAux<1,2> IPAux12;
    typedef IntersectionPointAux<1,3> IPAux13;
    typedef IntersectionPointAux<2,3> IPAux23;


    /** @brief Default constructor, creates empty object.
     * Resizes vectors for Plucker coordinates and products.
     */
    ComputeIntersection();

    /** @brief Constructor, sets both triangle objects.
     * It allocates memory, computes plucker coordinates and products.
     * @param triangle intersecting triangle object
     * @param tetrahedron intersecting tetrahedron object
     */
    ComputeIntersection(ElementAccessor<3> triangle, ElementAccessor<3> tetrahedron);
    ~ComputeIntersection();

    /** @brief Initializes lower dimensional objects.
     * Sets correctly the pointers to Plucker coordinates and products.
     */
    void init();

    /** @brief Computes intersection points for 2D-3D intersection.
     * Computes lower dimensional CIs:
     *      1) 3x triangle side vs tetrahedron
     *      2) 6x tetrahedron edge vs triangle
     * During cycle determine topology information and set computed flags accordingly.
     * Finally, tracing algorithm is called to sort the intersection points to
     * create convex polygon.
     * @param intersection final 2D-3D intersection object (output)
     *          It is then used further in prolongation decision routines.
     * @return number of intersection points found
     */
    void compute(IntersectionAux<2,3> &intersection);



    /// Prints out the Plucker coordinates of triangle sides and tetrahedron edges.
    void print_plucker_coordinates(std::ostream &os);
    /// Prints out the Plucker coordinates in a tree simplices.
    void print_plucker_coordinates_tree(std::ostream &os);

private:
    // S3 n-face pair (for edge-triangle phase)
    typedef std::array<uint, 2> FacePair;


    const unsigned int no_idx;
    // TODO rename s4 to s3, s3 to s2
    std::vector<unsigned int> s3_dim_starts; // get global index of n-face i of dimension d: s4_dim_starts[d]+i
    std::vector<unsigned int> s2_dim_starts;  // same for n-faces of S2

    std::vector<IPAux12> IP12s_;
    std::vector<IPAux23> IP23_list; //, degenerate_ips;

    /// Vector of Plucker coordinates for triangle side.
    std::vector<Plucker *> plucker_coordinates_triangle_;
    /// Vector of Plucker coordinates for tetrahedron edges.
    std::vector<Plucker *> plucker_coordinates_tetrahedron;

    /// Vector of pointers to Plucker products of triangle sides and tetrahedron edges.
    std::vector<double *> plucker_products_;
    
    /// Compute 1D-3D intersection objects [3]
    ComputeIntersection<1,3> CI13[3];
    /// Compute 1D-2D intersection objects [6]
    ComputeIntersection<1,2> CI12[6];


    // successors of IPs
    std::vector<unsigned int> IP_next;
    // successors of n-face objects
    // 4 vertices, 6 edges, 4 faces, 1 volume, 3 corners, 3 sides, 1 surface; total 22
    std::vector<unsigned int> object_next;

    //bool obj_have_back_link(unsigned int i_obj);
    auto edge_faces(uint i_edge) -> FacePair;

    auto vertex_faces(uint i_vtx) -> FacePair;

    /**
     * Returns true if: i_obj -> IP -> i_obj
     * Returns false if i_obj -> no_idx OR i_obj -> IP -> j_obj != i_obj
     *
     * Backlink marks: end of the created chain
     */
    inline bool have_backlink(uint i_obj);

    /**
     * Add raw intersection point.
     */
    inline unsigned int add_ip(const IPAux23 &ip) {
        //DebugOut() << "IP[" << IP23_list.size() << "]:" << ip;
        IP23_list.push_back(ip);
        return IP23_list.size() - 1;
    }

    /**
     * Set links: obj_before -> IP -> obj_after
     * if obj_after have null successor, set obj_after -> IP (backlink)
     */
    inline void set_links(uint obj_before_ip, uint ip_idx, uint obj_after_ip);
    
    const std::vector<std::vector<arma::uvec>> on_faces;
    std::vector<std::vector<arma::uvec>> _on_faces();

    bool S3_inverted;
    IntersectionAux<2,3>* intersection_;
};



#endif  // COMPUTE_INTERSECTION_H_
