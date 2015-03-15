/*
 * computeintersection.h
 *
 *
 *      Author: viktor
 */

//#include <armadillo>
#include "plucker.h"
#include "simplex.h"

//#include "intersectionpoint.h"
#include "intersectionlocal.h"
#include "system/system.hh"

using namespace std;
namespace computeintersection {


template<class A, class B> class ComputeIntersection {};

/******************************************************************
 * 	TŘÍDA PRO VÝPOČET SIMPLEX 1 - SIMPLEX 2
 * ****************************************************************/

template<> class ComputeIntersection<Simplex<1>, Simplex<2>> {


public:

	ComputeIntersection();
	ComputeIntersection(Simplex<1> &abs, Simplex<2> &triang);
	inline ~ComputeIntersection() {};

	void clear_all();
	bool compute(std::vector<IntersectionPoint<1,2>> &IP12s, bool compute_zeros_plucker_products);
	void initPluckerToCompute();
	void set_data(Simplex<1> *abs, Simplex<2> *triang);


	inline void set_pc_abscissa(Plucker *p){
		plucker_coordinates_abscissa[0] = p;
	}
	inline void set_pc_triangle(Plucker *p, unsigned int index){
		plucker_coordinates_triangle[index] = p;
	}

	inline Plucker *get_pc_abscissa(){
		return plucker_coordinates_abscissa[0];
	}

	inline Plucker *get_pc_triangle(unsigned int index){
		return plucker_coordinates_triangle[index];
	}

	void toStringPluckerCoordinates();

	void setPluckerProduct(double* number, unsigned int i);
	double* getPluckerProduct(unsigned int i);

	bool isComputed();
	void setComputed();


private:
	Simplex<1> *abscissa;
	Simplex<2> *triangle;

	std::vector<Plucker *> plucker_coordinates_abscissa;
	std::vector<Plucker *> plucker_coordinates_triangle;

	double *plucker_products[3];
	bool computed;

	static const double epsilon;

};

/******************************************************************
 * 	TŘÍDA PRO VÝPOČET SIMPLEX 1 - SIMPLEX 3
 * ****************************************************************/

template<> class ComputeIntersection<Simplex<1>, Simplex<3>> {

public:

	ComputeIntersection();
	ComputeIntersection(Simplex<1> &abs,Simplex<3> &tetr);

	void clear_all();
	void init();
	void set_data(Simplex<1> *abs, Simplex<3> *tetr);
	int compute(std::vector<IntersectionPoint<1,3>> &IP13s);

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

	void toStringPluckerCoordinates();
	void toStringPluckerCoordinatesTree();

	void setPluckerProduct(double* number, unsigned int index_CI, unsigned index_edge);
	double* getPluckerProduct(unsigned int index_CI, unsigned index_edge);

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

	ComputeIntersection(Simplex<2> &triangle, Simplex<3> &tetr);

	void clear_all();
	void init();
	void compute(IntersectionLocal &lokalni_mnohouhlenik);

	void toStringPluckerCoordinates();
	void toStringPluckerCoordinatesTree();

	inline ~ComputeIntersection() {};

//private:

	// Reprezentation of triangle and tetrahedron as object Simplex
	Simplex<2> *triange;
	Simplex<3> *tetrahedron;

	// Plucker coordinates for each abscissa of simplices
	std::vector<Plucker *> plucker_coordinates_triangle;
	std::vector<Plucker *> plucker_coordinates_tetrahedron;

	// Computing objects
	ComputeIntersection<Simplex<1>, Simplex<3>> CI13[3];
	ComputeIntersection<Simplex<1>, Simplex<2>> CI12[6];

	// Intersection objects
	//std::vector<IntersectionLocal> intersections;


};

} // END namespace_close
