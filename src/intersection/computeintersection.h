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

namespace computeintersection {

// forward declare
template<class A, class B> class ComputeIntersection;
class IntersectionPolygon;
class Plucker;
template<unsigned int, unsigned int> class IntersectionPoint;

/******************************************************************
 * 	TŘÍDA PRO VÝPOČET SIMPLEX 1 - SIMPLEX 2
 * ****************************************************************/

/**
 * TODO: comment
 */
template<> class ComputeIntersection<Simplex<1>, Simplex<2>> {
public:

	ComputeIntersection();
    //TODO: pass object by pointers or const reference, unify with set_data()
	ComputeIntersection(Simplex<1> &abs, Simplex<2> &triang);
	inline ~ComputeIntersection() {};

	
    //TODO: why this is not done in constructor?
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
	bool compute(std::vector<IntersectionPoint<1,2>> &IP12s, bool compute_zeros_plucker_products);
    
    
	void set_data(Simplex<1> *abs, Simplex<2> *triang);


    /// Sets the pointer to Plucker coordinates of the abscissa.
	inline void set_pc_abscissa(Plucker *p){
		plucker_coordinates_abscissa[0] = p;
	}
	
	/// Sets the pointer to Plucker coordinates of the triangle side of given @p side_idx.
	inline void set_pc_triangle(Plucker *p, unsigned int side_idx){
		plucker_coordinates_triangle[side_idx] = p;
	}

	/// Gets the pointer to Plucker coordinates of the abscissa.
	inline Plucker *get_pc_abscissa(){
		return plucker_coordinates_abscissa[0];
	}

	/// Gets the pointer to Plucker coordinates of the triangle side of given @p side_idx.
	inline Plucker *get_pc_triangle(unsigned int side_idx){
		return plucker_coordinates_triangle[side_idx];
	}

	void print_plucker_coordinates(std::ostream &os);

	void set_plucker_product(double* number, unsigned int i);
	double* get_plucker_product(unsigned int i);

	bool is_computed();
	void set_computed();


private:
    /// Resets Plucker products to NULL.
    void clear_all();
    /// Computes Plucker coordinates (abscissa, triangle lines) and Plucker products.
    void compute_plucker_products();
    
	Simplex<1> *abscissa;
	Simplex<2> *triangle;

    //TODO: Is there a reason for abscissa pc to be a vector?
	std::vector<Plucker *> plucker_coordinates_abscissa;
	std::vector<Plucker *> plucker_coordinates_triangle;

    //TODO: allocate at the top level intersection object, use NaN to indicate plucker products not computed yet, also Destroy!
	double *plucker_products[3];
	bool computed;
};

/******************************************************************
 * 	TŘÍDA PRO VÝPOČET SIMPLEX 1 - SIMPLEX 3
 * ****************************************************************/

template<> class ComputeIntersection<Simplex<1>, Simplex<3>> {

public:

	ComputeIntersection();
	ComputeIntersection(Simplex<1> &abs,Simplex<3> &tetr);

	
	void init();
	void set_data(Simplex<1> *abs, Simplex<3> *tetr);
    //TODO comment cases in implementation
	unsigned int compute(std::vector<IntersectionPoint<1,3>> &IP13s);

	inline void set_pc_abscissa(Plucker *p){
		plucker_coordinates_abscissa[0] = p;
	}
	inline void set_pc_tetrahedron(Plucker *p, unsigned int index){
		plucker_coordinates_tetrahedron[index] = p;
	}

	inline Plucker *get_pc_abscissa(){
		return plucker_coordinates_abscissa[0];
	}

	inline Plucker *get_pc_tetrahedron(unsigned int index){
		return plucker_coordinates_tetrahedron[index];
	}

	void print_plucker_coordinates(std::ostream &os);
	void print_plucker_coordinates_tree(std::ostream &os);

	void set_plucker_product(double* number, unsigned int index_CI, unsigned index_edge);
	double* get_plucker_product(unsigned int index_CI, unsigned index_edge);

	inline ~ComputeIntersection() {}

private:
	Simplex<1> *abscissa;
	Simplex<3> *tetrahedron;

	std::vector<Plucker *> plucker_coordinates_abscissa;
	std::vector<Plucker *> plucker_coordinates_tetrahedron;

	ComputeIntersection<Simplex<1>, Simplex<2>> CI12[4];
};

/******************************************************************
 * 	TŘÍDA PRO VÝPOČET SIMPLEX 2 - SIMPLEX 3
 * ****************************************************************/

template<> class ComputeIntersection<Simplex<2>, Simplex<3> > {
public:
	ComputeIntersection();

	ComputeIntersection(Simplex<2> &tria, Simplex<3> &tetr);

	void init();
	void compute(IntersectionPolygon &lokalni_mnohouhlenik);

	void print_plucker_coordinates(std::ostream &os);
	void print_plucker_coordinates_tree(std::ostream &os);

	inline ~ComputeIntersection() {};

private:

	// Representation of triangle and tetrahedron as object Simplex
	Simplex<2> *triangle;
	Simplex<3> *tetrahedron;

	// Plucker coordinates for each abscissa of simplices
	std::vector<Plucker *> plucker_coordinates_triangle;
	std::vector<Plucker *> plucker_coordinates_tetrahedron;

	// Computing objects
	ComputeIntersection<Simplex<1>, Simplex<3>> CI13[3];
	ComputeIntersection<Simplex<1>, Simplex<2>> CI12[6];
};

} // END namespace_close


#endif  // COMPUTE_INTERSECTION_H_