/*
 *      Author: viktor
 */

#include "computeintersection.h"
#include "system/system.hh"
namespace computeintersection{

Plucker getPluckerFromSimplex(const Simplex<1> &abs){
	Plucker pl(abs[0].getPointCoordinates(),abs[1].getPointCoordinates());
	return pl;
};

Simplex<1> getAbsicca(const Simplex<2> &abs, int i){
	return abs[i];
};

Simplex<1> getAbsicca(const Simplex<3> &abs, int i){
	if(i < 3){
		return abs[0][i];
	}else if(i < 5){
		return abs[1][i - 2];
	}else{
		return abs[2][2];
	}
};


/**********************************************
 * METODY PRO SIMPLEX 1 A SIMPLEX 3
 * ********************************************/
/*
ComputeIntersection<Simplex<1>, Simplex<3>>::ComputeIntersection(){
	xprintf(Msg, "Konstruktor - prazdny - 13\n");
};

ComputeIntersection<Simplex<1>, Simplex<3>>::ComputeIntersection(const Simplex<1> &abs, const Simplex<3> &tetrahedron){
	this->abscissa = abs;
	this->tetrahedron = tetrahedron;
	//xprintf(Msg, "Konstruktor - s parametry 13\n");
};

void ComputeIntersection<Simplex<1>, Simplex<3>>::compute(){
	int smer_minus = -1;
	int smer_plus = 1;
*//*
	//Optimalizovaná verze
	ComputeIntersection<Simplex<1>, Simplex<2>> sim_0(abscissa, tetrahedron[0]);
	sim_0.setPluckerDirection(smer_minus,1);
	sim_0.compute();
	ComputeIntersection<Simplex<1>, Simplex<2>> sim_1(abscissa, tetrahedron[1]);
	sim_1.setPluckerDirection(smer_minus,0);
	sim_1.setPluckerDirection(smer_minus,2);
	sim_1.setPluckerProduct(sim_0.getPluckerProduct(0),0);
	sim_1.compute();
	ComputeIntersection<Simplex<1>, Simplex<2>> sim_2(abscissa, tetrahedron[2]);
	sim_2.setPluckerDirection(smer_minus,1);
	sim_2.setPluckerProduct(sim_0.getPluckerProduct(1),0);
	sim_2.setPluckerProduct(sim_1.getPluckerProduct(1),1);
	sim_2.compute();
	ComputeIntersection<Simplex<1>, Simplex<2>> sim_3(abscissa, tetrahedron[3]);
	sim_3.setPluckerDirection(smer_minus,0);
	sim_3.setPluckerDirection(smer_minus,2);
	sim_3.setPluckerProduct(sim_0.getPluckerProduct(2),0);
	sim_3.setPluckerProduct(sim_1.getPluckerProduct(2),1);
	sim_3.setPluckerProduct(sim_2.getPluckerProduct(2),2);
	sim_3.compute();*/

	// Neoptimalizovana verze:
	/*for(int i = 0; i < 4; i++){
		ComputeIntersection<Simplex<1>, Simplex<2>> sim12(abscissa, tetrahedron[i]);
		sim12.compute();
	}*/

	//xprintf(Msg, "Hell yeah - spousti sa \n");
//};

/****************************************************************
 * METODY PRO SIMPLEX 1 A SIMPLEX 2
 ****************************************************************/
/*
ComputeIntersection<Simplex<1>, Simplex<2>>::ComputeIntersection(){
	intersectionexists = false;
};

ComputeIntersection<Simplex<1>, Simplex<2>>::ComputeIntersection(const Simplex<1> &abs, const Simplex<2> &triang){
	intersectionexists = false;
	abscissa = abs;
	triangle = triang;
	plucker_products[0] = NULL;
	plucker_products[1] = NULL;
	plucker_products[2] = NULL;
	plucker_directions[0] = 1;
	plucker_directions[1] = 1;
	plucker_directions[2] = 1;
};

void ComputeIntersection<Simplex<1>, Simplex<2>>::setSimplices(const Simplex<1> &abs, const Simplex<2> &triang){
	abscissa = abs;
	triangle = triang;
};
void ComputeIntersection<Simplex<1>, Simplex<2>>::setPluckerProduct(double &number, int i){
	plucker_products[i] = &number;
};
double &ComputeIntersection<Simplex<1>, Simplex<2>>::getPluckerProduct(int i){
	return *plucker_products[i];
}

void ComputeIntersection<Simplex<1>, Simplex<2>>::setPluckerDirection(int &number, unsigned int i){
	plucker_directions[i] = number;
};
int &ComputeIntersection<Simplex<1>, Simplex<2>>::getPluckerDirection(unsigned int i){
	return plucker_directions[i];
};

Simplex<1> &ComputeIntersection<Simplex<1>, Simplex<2>>::getFirstSimplex(){
	return abscissa;
};
Simplex<2> &ComputeIntersection<Simplex<1>, Simplex<2>>::getSecondSimplex(){
	return triangle;
};

bool ComputeIntersection<Simplex<1>, Simplex<2>>::intersection_exists(){
	return intersectionexists;
};

bool ComputeIntersection<Simplex<1>, Simplex<2>>::getDirection(){
	return abscissa_direction;
};

void ComputeIntersection<Simplex<1>, Simplex<2>>::compute(){

	for(int i = 0; i < 3; i++){
		if(plucker_products[i] == NULL){
		plucker_products[i] = new double(getPluckerFromSimplex(abscissa)*getPluckerFromSimplex(getAbsicca(triangle,i)));
		cout << "PP nastaven = " << *plucker_products[i] << " smer: "<<  plucker_directions[i] << endl;
		}
	}

	// Řešit orientace, REFF elementy atd.

	if((*plucker_products[0]*plucker_directions[0] > 0) && (*plucker_products[1]*plucker_directions[1] > 0) && (*plucker_products[2]*plucker_directions[2] > 0)){
		cout << "Hell Yeah Průnik" << endl;
		// ukládání Intersection local -> nebo jeho vracení radši
	}else if((*plucker_products[0]*plucker_directions[0] < 0) && (*plucker_products[1]*plucker_directions[1] < 0) && (*plucker_products[2]*plucker_directions[2] < 0)){
		cout << "Hell Yeah Průnik" << endl;
		// ukládání Intersection local -> nebo jeho vracení radši
	}else{
		// Vracet NULL obejkt
	}
};*/
/* KONEC metoda pro simplex 1 a simplex 2 */


/****************************************************************
 * METODY PRO SIMPLEX 1 A SIMPLEX 2
 ****************************************************************/
ComputeIntersection<Simplex<1>, Simplex<2>>::ComputeIntersection(){};


/****************************************************************
 * METODY PRO SIMPLEX 1 A SIMPLEX 3
 ****************************************************************/
ComputeIntersection<Simplex<1>, Simplex<3>>::ComputeIntersection(){};

ComputeIntersection<Simplex<1>, Simplex<3>>::ComputeIntersection(const Simplex<1> &abscissa, const Simplex<3> &tetrahedron){

};


/****************************************************************
 * METODY PRO SIMPLEX 2 A SIMPLEX 3
 ****************************************************************/
ComputeIntersection<Simplex<2>, Simplex<3>>::ComputeIntersection(){};


ComputeIntersection<Simplex<2>, Simplex<3>>::ComputeIntersection(const Simplex<2> &triangle, const Simplex<3> &tetrahedron){
	this->triange = triangle;
	this->tetrahedron = tetrahedron;

	this->clear_all();

	for(unsigned int i = 0; i < 6;i++){
		//this->CI12[i] = ComputeIntersection<Simplex<1>, Simplex<2>>(tetrahedron.getAbscissa(i) ,triangle);
	}
	for(unsigned int i = 0; i < 3;i++){
		this->CI13[i] = ComputeIntersection<Simplex<1>, Simplex<3>>(triangle[i] ,tetrahedron);
		// nebo když bude CI13 pointer
		// this->CI13[i] = new ComputeIntersection<Simplex<1>, Simplex<3>>(triangle[i] ,tetrahedron);
	}

};

void ComputeIntersection<Simplex<2>, Simplex<3>>::clear_all(){
	for(unsigned int i = 0; i < 6;i++){
		this->p_coordinates_tetrahedron[i] = NULL;
		//this->CI12[i] = NULL;
	}
	for(unsigned int i = 0; i < 3;i++){
		this->p_coordinates_triange[i] = NULL;
		//this->CI13[i] = NULL;
	}
};

void ComputeIntersection<Simplex<2>, Simplex<3>>::init(){

};

void ComputeIntersection<Simplex<2>, Simplex<3>>::compute(){

};

void ComputeIntersection<Simplex<2>, Simplex<3>>::setPC_triangle(Plucker* p_triangle_coordinate, unsigned int index){
	p_coordinates_triange[index] = p_triangle_coordinate;
};

void ComputeIntersection<Simplex<2>, Simplex<3>>::setPC_tetrahedron(Plucker* p_tetrahedron_coordinate, unsigned int index){
	p_coordinates_tetrahedron[index] = p_tetrahedron_coordinate;
};

} // END namespace_close
