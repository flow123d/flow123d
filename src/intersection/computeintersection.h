/*
 * computeintersection.h
 *
 *
 *      Author: viktor
 */

//#include <armadillo>
#include "plucker.h"
#include "simplex.h"
#include "system/system.hh"

using namespace std;
namespace computeintersection {
// Metody pro konvezi objektů

Plucker getPluckerFromSimplex(const Simplex<1> &abs);
Simplex<1> getAbsicca(const Simplex<2> &abs, int i);
Simplex<1> getAbsicca(const Simplex<3> &abs, int i);

// Výpočetní třídy

template<class A, class B> class ComputeIntersection {};

/******************************************************************
 * 	TŘÍDA PRO VÝPOČET SIMPLEX 1 - SIMPLEX 2
 * ****************************************************************/

template<> class ComputeIntersection<Simplex<1>, Simplex<2>> {


public:

	ComputeIntersection();
	ComputeIntersection(const Simplex<1> &abs, const Simplex<2> &triang);
	inline ~ComputeIntersection() {}

	void compute();

	void setPluckerProduct(double &number, int i);
	double &getPluckerProduct(int i);

	//void setPluckerDirection(int &number, unsigned int i);
	//int &getPluckerDirection(unsigned int i);

	//bool intersection_exists();
	//bool getDirection();

private:
	Simplex<1> abscissa;
	Simplex<2> triangle;

	// Pro tento typ průniku existují právě 3 plückerovy součiny
	double* plucker_products[3];

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
	ComputeIntersection(const Simplex<1> &abscissa, const Simplex<3> &tetrahedron);

	void init();
	void compute();

	void setPC_abscissa(const Plucker &p_abscissa_coordinate);
	void setPC_tetrahedron(const Plucker &p_tetrahedron_coordinate, unsigned int index);

	inline ~ComputeIntersection() {}

private:
	Simplex<1> abscissa;
	Simplex<3> tetrahedron;

	Plucker *p_coordinates_abscissa[1];
	Plucker *p_coordinates_tetrahedron[6];

	ComputeIntersection<Simplex<1>, Simplex<2>> CI12[4];

};

/******************************************************************
 * 	TŘÍDA PRO VÝPOČET SIMPLEX 2 - SIMPLEX 3
 * ****************************************************************/

template<> class ComputeIntersection<Simplex<2>, Simplex<3> > {
public:
	ComputeIntersection();

	ComputeIntersection(const Simplex<2> &triangle,
						const Simplex<3> &tetrahedron);

	void clear_all();
	void init();
	void compute();

	void setPC_triangle(Plucker* p_triangle_coordinate, unsigned int index);
	void setPC_tetrahedron(Plucker* p_tetrahedron_coordinate, unsigned int index);

	inline ~ComputeIntersection() {}

//private:

	// Reprezentation of triangle and tetrahedron as object Simplex
	Simplex<2> triange;
	Simplex<3> tetrahedron;

	// Plucker coordinates for each abscissa of simplices
	/* pozn. Myslet na to, že bude jakýsi globalní Plucker, který si bude načítat už vytvořená data
	 * a vracet vždy stejný objekt, jen s jinými daty*/
	Plucker *p_coordinates_triange[3];
	Plucker *p_coordinates_tetrahedron[6];

	// Computing objects
	ComputeIntersection<Simplex<1>, Simplex<3>> CI13[3];
	ComputeIntersection<Simplex<1>, Simplex<2>> CI12[6];

};

} // END namespace_close
