/*
 * computeintersection.h
 *
 *
 *      Author: viktor
 */
#ifndef COMPUTE_INTERSECTION_H_
#define COMPUTE_INTERSECTION_H_

#include "simplex.h"
#include "system/system.hh"
#include "mesh/ref_element.hh"

namespace computeintersection {

// forward declare
template<class A, class B> class ComputeIntersection;
class Plucker;
template<unsigned int, unsigned int> class IntersectionAux;
template<unsigned int, unsigned int> class IntersectionPointAux;

static const double plucker_empty = std::numeric_limits<double>::infinity();

/******************************************************************
 * 	TŘÍDA PRO VÝPOČET SIMPLEX 1 - SIMPLEX 2
 * ****************************************************************/

/**
 * TODO: comment
 * TODO: Idea: constructor creates empty object; for each abscissa and triangle only update data 
 * [abscissa, triangle] and (plucker coordinates, products) if not final dimension.
 * Remove init() methods.
 * Wouldn't a flag for final dimension object help code structure?
 * 
 * When solve allocating plucker products, solve also destroying plucker coordinates.
 * 
 * @brief Computes the intersection of an abscissa and a triangle.
 */
template<> class ComputeIntersection<Simplex<1>, Simplex<2>> {
public:
    /// Default constructor. Use when this is NOT final intersection object.
	ComputeIntersection();
    /** @brief Constructor, sets abscissa and triangle object.
     * Use when this is final intersection object.
     * It allocates memory, computes plucker coordinates and products.
     */
	ComputeIntersection(Simplex<1> &abscissa, Simplex<2> &triangle);
	~ComputeIntersection();

    // Why this is not done in constructor?
    // Because default constructor is called in 1d-3d, 2d-3d and compute() is called later.
    
    /** @brief Computes intersection points of line and triangle.
     * 
     * If Plucker products are nonezero and with the same sign, then IP is inside the triangle.
     * If some of the Plucker products are zero:
     * 1 zero product -> IP is on the triangle side
     * 2 zero products -> IP is at the vertex of triangle (there is no other IP)
     * 3 zero products: 
     *      -> IP is at the vertex of triangle but the line is parallel to opossite triangle side
     *      -> triangle side is part of the line (and otherwise)     
     * IP is intersection of triangle and whole line (bisector).
     * @param IP12s - input/output vector of IPs. If IP found, it is pushed back.
     * @param compute_zeros_plucker_products - if true, resolve pathologic cases (zero Plucker products), 
     * otherwise ignore. E.g. in 2d-3d is false when looking for tetrahedron edges X triangle intersection
     * (these would be found before in triangle line X tetrahedron intersection).
     * @return true, if intersection is found; false otherwise
     */
	bool compute(std::vector<IntersectionPointAux<1,2>> &IP12s, bool compute_zeros_plucker_products);
    
    /** Computes final 1d-2d intersection. (Use when this is not the resulting dimension object).
     * @param IP12s - input/output vector of IPs. If IP found, it is pushed back.
     * @return true, if intersection is found; false otherwise
     */
    bool compute_final(std::vector<IntersectionPointAux<1,2>> &IP12s);
    
    /// @name Setters and Getters
    //@{ 
    /** @brief Sets the abscissa and triangle.
     * 
     * Use mainly when this is not final intersection computation.
     */
	void set_data(Simplex<1> &abscissa, Simplex<2> &triangle);

    /// Sets the pointer to Plucker coordinates of the abscissa.
    void set_pc_abscissa(Plucker *p){
		plucker_coordinates_abscissa_ = p;
	}
	
	/// Sets the pointer to Plucker coordinates of the triangle side of given @p side_idx.
    void set_pc_triangle(Plucker *p, unsigned int side_idx){
		plucker_coordinates_triangle_[side_idx] = p;
	}

	/// Gets the pointer to Plucker coordinates of the abscissa.
    Plucker *get_pc_abscissa(){
		return plucker_coordinates_abscissa_;
	}

	/// Gets the pointer to Plucker coordinates of the triangle side of given @p side_idx.
    Plucker *get_pc_triangle(unsigned int side_idx){
		return plucker_coordinates_triangle_[side_idx];
	}

	/** @brief Sets the Plucker product of abscissa and triangle side (use if computed before).
     * @param side_idx local index of a side of triangle.
     */
    void set_plucker_product(double* number, unsigned int side_idx)
    { plucker_products_[side_idx] = number; };
    
    /** @brief Getter for Plucker product of abscissa and triangle side @p side_idx.
     */
    double* get_plucker_product(unsigned int side_idx)
    { return plucker_products_[side_idx]; };

    /// Gets true if the intersection has been computed already (e.g. in case of IP in vertex).
    bool is_computed() { return computed_; }
    /// Sets the 'computed' flag true. Means that intersection has been computed already (e.g. in case of IP in vertex).
    void set_computed() { computed_ = true; }
    
	//@}

    /// Prints out all the Plucker coordinates.
	void print_plucker_coordinates(std::ostream &os);

    /// Resets Plucker products to 'nullptr'.
    /// Use this CAREFULLY as it does not destroy the objects.
    /// It is intended to be used only from higher dimensions when before destroying.
    void clear_all();
    
private:
    /// Computes Plucker coordinates (abscissa, triangle lines) and Plucker products.
    void compute_plucker_products();
    
    double signed_plucker_product(unsigned int i)
    { return RefElement<2>::normal_orientation(i) ? -(*plucker_products_[i]) : (*plucker_products_[i]); }
    
        /** Computes intersection when nonezero Plucker products are of the same sign.
     * @param IP is the intersection point (if found)
     * @return true, if intersection is found; false otherwise
     */
    bool compute_plucker(IntersectionPointAux<1,2> &IP);
    
    /** Computes intersection of abscissa and triangle side for zero Plucker product - pathologic case.
     * @param side is the local index of the triangle side
     * @param IP is the intersection point (if found)
     * @return true, if intersection is found; false otherwise
     */
    bool compute_pathologic(unsigned int side, IntersectionPointAux<1,2> &IP);
    
    /// Flag 'computed'; is true is intersection has been computed already.
    bool computed_;
    
	Simplex<1> *abscissa_;
	Simplex<2> *triangle_;

	Plucker* plucker_coordinates_abscissa_;
	std::vector<Plucker *> plucker_coordinates_triangle_;
    /// Pointers to Plucker products of abscissa and triangle side.
	std::vector<double *> plucker_products_;
};



/******************************************************************
 * 	TŘÍDA PRO VÝPOČET SIMPLEX 1 - SIMPLEX 3
 * ****************************************************************/

/**
 * @brief Computes the intersection of an abscissa and a tetrahedron.
 * 
 * Uses 4 ComputeIntersection<Simplex<1>,Simplex<2>> for abscissa and tetrahedron sides intersections.
 */
template<> class ComputeIntersection<Simplex<1>, Simplex<3>> {

public:

	/// Default constructor. Use when this is NOT final intersection object.
    ComputeIntersection();
    /** @brief Constructor, sets abscissa and tetrahedron object.
     * Use when this is final intersection object.
     * It allocates memory, computes plucker coordinates and products.
     */
	ComputeIntersection(Simplex<1> &abscissa,Simplex<3> &tetrahedron);
	~ComputeIntersection();
	
	void init();
    
    //TODO comment cases in implementation
    unsigned int compute(std::vector<IntersectionPointAux<1,3>> &IP13s);
    unsigned int compute(IntersectionAux<1,3> &intersection, std::vector<unsigned int> &prolongation_table);
    
     /// @name Setters and Getters
    //@{ 
    /**
     * @brief Sets the abscissa and tetrahedron.
     * Use mainly when this is not final intersection computation.
     */
    void set_data(Simplex<1> &abscissa, Simplex<3> &tetrahedron);
    
    /// Sets the pointer to Plucker coordinates of the abscissa.
    void set_pc_abscissa(Plucker *p){
		plucker_coordinates_abscissa_ = p;
	}
	
	/// Sets the pointer to Plucker coordinates of the tetrahedron edge of given @p edge_idx.
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
    
    /// Sets the pointer to Plucker product of abscissa and tetrahedron edge of given @p edge_idx.
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
    
    Plucker* plucker_coordinates_abscissa_;
	std::vector<Plucker *> plucker_coordinates_tetrahedron;
    /// Pointers to Plucker products of abscissa and tetrahedron edges.
    std::vector<double *> plucker_products_;
	ComputeIntersection<Simplex<1>, Simplex<2>> CI12[4];
};

/******************************************************************
 * 	TŘÍDA PRO VÝPOČET SIMPLEX 2 - SIMPLEX 3
 * ****************************************************************/

/**
 * @brief Computes the intersection of a triangle and a tetrahedron.
 * 
 * Uses 3 ComputeIntersection<Simplex<1>,Simplex<3>> for triangle sides and tetrahedron intersections.
 * Uses 6 ComputeIntersection<Simplex<1>,Simplex<2>> for tetrahedron sides and triangle intersections.
 */
template<> class ComputeIntersection<Simplex<2>, Simplex<3> > {
public:
	ComputeIntersection();

	ComputeIntersection(Simplex<2> &triangle, Simplex<3> &tetrahedron);
    ~ComputeIntersection();

	void init();
// 	void compute(IntersectionPolygon &local_polygon);
//    void compute(std::vector<IntersectionPointAux<2,3>> &IP23s);
    /// @brief 
    /** 
     * @param prolongation_table is an auxiliary vector that is filled in tracing algorithm of polygon.
     * It is then used further in prolongation decision routines.
     */
    void compute(IntersectionAux<2,3> &intersection, std::vector<unsigned int> &prolongation_table);

    /// Prints out the Plucker coordinates of triangle sides and tetrahedron edges.
	void print_plucker_coordinates(std::ostream &os);
    /// Prints out the Plucker coordinates in a tree simplices.
	void print_plucker_coordinates_tree(std::ostream &os);

private:

	// Plucker coordinates for each abscissa of simplices
	std::vector<Plucker *> plucker_coordinates_triangle_;
	std::vector<Plucker *> plucker_coordinates_tetrahedron;

    /// Pointers to Plucker products of triangle sides and tetrahedron edges.
    std::vector<double *> plucker_products_;
    
	// Computing objects
	ComputeIntersection<Simplex<1>, Simplex<3>> CI13[3];
	ComputeIntersection<Simplex<1>, Simplex<2>> CI12[6];
};

} // END namespace_close


#endif  // COMPUTE_INTERSECTION_H_