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
// Metody pro konvezi objektů

//Plucker getPluckerFromSimplex(const Simplex<1> &abs);
//Simplex<1> getAbsicca(const Simplex<2> &abs, int i);
//Simplex<1> getAbsicca(const Simplex<3> &abs, int i);

// Výpočetní třídy

template<class A, class B> class ComputeIntersection {};

/******************************************************************
 * 	TŘÍDA PRO VÝPOČET SIMPLEX 1 - SIMPLEX 2
 * ****************************************************************/

template<> class ComputeIntersection<Simplex<1>, Simplex<2>> {


public:

	ComputeIntersection();
	ComputeIntersection(Simplex<1> &abs, Simplex<2> &triang);
	//ComputeIntersection(Simplex<1> &abs,Simplex<2> &triang);
	inline ~ComputeIntersection() {};

	void clear_all();
	bool compute(IntersectionPoint<1,2> &IP, bool compute_zeros_plucker_products);
	void initPluckerToCompute();


	std::vector<Plucker *> &getPC_abscissa();
	std::vector<Plucker *> &getPC_triangle();
	Plucker &getPC_abscissa(unsigned int index);
	Plucker &getPC_triangle(unsigned int index);

	void setPC_abscissa(std::vector<Plucker *> &p_abscissa_coordinates);
	void setPC_triangle(std::vector<Plucker *> &p_triangle_coordinates);
	void setPC_abscissa(Plucker &p_abscissa_coordinate);
	void setPC_triangle(Plucker &p_triangle_coordinate, unsigned int index);

	void toStringPluckerCoordinates();

	void setPluckerProduct(double* number, unsigned int i);
	double* getPluckerProduct(unsigned int i);

	bool isComputed();
	void setComputed();
	//void setPluckerDirection(int &number, unsigned int i);
	//int &getPluckerDirection(unsigned int i);

	//bool intersection_exists();
	//bool getDirection();

private:
	Simplex<1> *abscissa;
	Simplex<2> *triangle;

	std::vector<Plucker *> plucker_coordinates_abscissa;
	std::vector<Plucker *> plucker_coordinates_triangle;

	double *plucker_products[3];
	bool computed;


	//int plucker_directions[3];
	//bool intersectionexists;
	//int abscissa_direction;
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
	int compute(std::vector<IntersectionPoint<1,3>> &IP13s);

	std::vector<Plucker *> &getPC_abscissa();
	std::vector<Plucker *> &getPC_tetrahedron();
	Plucker &getPC_abscissa(unsigned int index);
	Plucker &getPC_tetrahedron(unsigned int index);

	void setPC_abscissa(std::vector<Plucker *> &p_abscissa_coordinates);
	void setPC_tetrahedron(std::vector<Plucker *> &p_tetrahedron_coordinates);
	void setPC_abscissa(Plucker &p_abscissa_coordinate);
	void setPC_tetrahedron(Plucker &p_tetrahedron_coordinate, unsigned int index);

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


	std::vector<Plucker *> &getPC_triangle();
	std::vector<Plucker *> &getPC_tetrahedron();
	Plucker &getPC_triangle(unsigned int index);
	Plucker &getPC_tetrahedron(unsigned int index);

	void setPC_triangle(std::vector<Plucker *> &p_triangle_coordinates);
	void setPC_tetrahedron(std::vector<Plucker *> &p_tetrahedron_coordinates);
	void setPC_triangle(Plucker &p_triangle_coordinate, unsigned int index);
	void setPC_tetrahedron(Plucker &p_tetrahedron_coordinate, unsigned int index);

	void toStringPluckerCoordinates();
	void toStringPluckerCoordinatesTree();

	inline ~ComputeIntersection() {}

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
