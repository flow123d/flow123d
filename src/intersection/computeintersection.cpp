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

ComputeIntersection<Simplex<1>, Simplex<2>>::ComputeIntersection(const Simplex<1> &abs,const Simplex<2> &triang){
	this->abscissa = abs;
	this->triangle = triang;

	this->clear_all();

};

void ComputeIntersection<Simplex<1>, Simplex<2>>::clear_all(){
	for(unsigned int i = 0; i < 3;i++){
		this->p_coordinates_triangle[i] = NULL;
	}
	this->p_coordinates_abscissa[0] = NULL;

};

void ComputeIntersection<Simplex<1>, Simplex<2>>::compute(){

};

Plucker* ComputeIntersection<Simplex<1>, Simplex<2>>::getPC_abscissa(){
	return p_coordinates_abscissa[0];
};

Plucker* ComputeIntersection<Simplex<1>, Simplex<2>>::getPC_triangle(unsigned int index){
	return p_coordinates_triangle[index];
};

void ComputeIntersection<Simplex<1>, Simplex<2>>::setPC_abscissa(Plucker* p_abscissa_coordinate){
	p_coordinates_abscissa[0] = p_abscissa_coordinate;
};

void ComputeIntersection<Simplex<1>, Simplex<2>>::setPC_triangle(Plucker* p_triangle_coordinate, unsigned int index){
	p_coordinates_triangle[index] = p_triangle_coordinate;
};



/****************************************************************
 * METODY PRO SIMPLEX 1 A SIMPLEX 3
 ****************************************************************/
ComputeIntersection<Simplex<1>, Simplex<3>>::ComputeIntersection(){};

ComputeIntersection<Simplex<1>, Simplex<3>>::ComputeIntersection(const Simplex<1> &abscissa,const Simplex<3> &tetrahedron){
	this->abscissa = abscissa;
	this->tetrahedron = tetrahedron;

	this->clear_all();

	for(unsigned int i = 0; i < 4;i++){
	CI12[i] = ComputeIntersection<Simplex<1>, Simplex<2>>(abscissa ,tetrahedron[i]);
	}
};

void ComputeIntersection<Simplex<1>, Simplex<3>>::clear_all(){
	for(unsigned int i = 0; i < 6;i++){
		this->p_coordinates_tetrahedron[i] = NULL;
	}
	this->p_coordinates_abscissa[0] = NULL;

};

void ComputeIntersection<Simplex<1>, Simplex<3>>::init(){
	for(unsigned int i = 0; i < 4;i++){
			CI12[i].setPC_abscissa(this->p_coordinates_abscissa[i]);
	 }
	CI12[0].setPC_triangle(this->p_coordinates_tetrahedron[0], 0);
	CI12[0].setPC_triangle(this->p_coordinates_tetrahedron[1], 1);
	CI12[0].setPC_triangle(this->p_coordinates_tetrahedron[2], 2);
	CI12[1].setPC_triangle(this->p_coordinates_tetrahedron[0], 0);
	CI12[1].setPC_triangle(this->p_coordinates_tetrahedron[3], 1);
	CI12[1].setPC_triangle(this->p_coordinates_tetrahedron[4], 2);
	CI12[2].setPC_triangle(this->p_coordinates_tetrahedron[1], 0);
	CI12[2].setPC_triangle(this->p_coordinates_tetrahedron[3], 1);
	CI12[2].setPC_triangle(this->p_coordinates_tetrahedron[5], 2);
	CI12[3].setPC_triangle(this->p_coordinates_tetrahedron[2], 0);
	CI12[3].setPC_triangle(this->p_coordinates_tetrahedron[4], 1);
	CI12[3].setPC_triangle(this->p_coordinates_tetrahedron[5], 2);

};

void ComputeIntersection<Simplex<1>, Simplex<3>>::compute(){

};

Plucker* ComputeIntersection<Simplex<1>, Simplex<3>>::getPC_abscissa(){
	return p_coordinates_abscissa[0];
};

Plucker* ComputeIntersection<Simplex<1>, Simplex<3>>::getPC_tetrahedron(unsigned int index){
	return p_coordinates_tetrahedron[index];
};

void ComputeIntersection<Simplex<1>, Simplex<3>>::setPC_abscissa(Plucker* p_abscissa_coordinate){
	p_coordinates_abscissa[0] = p_abscissa_coordinate;
};

void ComputeIntersection<Simplex<1>, Simplex<3>>::setPC_tetrahedron(Plucker* p_tetrahedron_coordinate, unsigned int index){
	p_coordinates_tetrahedron[index] = p_tetrahedron_coordinate;
};

/****************************************************************
 * METODY PRO SIMPLEX 2 A SIMPLEX 3
 ****************************************************************/
ComputeIntersection<Simplex<2>, Simplex<3>>::ComputeIntersection(){};


ComputeIntersection<Simplex<2>, Simplex<3>>::ComputeIntersection(Simplex<2> &triangle, Simplex<3> &tetrahedron){
	this->triange = triangle;
	this->tetrahedron = tetrahedron;

	this->clear_all();

	for(unsigned int i = 0; i < 6;i++){
		CI12[i] = ComputeIntersection<Simplex<1>, Simplex<2>>(tetrahedron.getAbscissa(i) ,triangle);
	}
	for(unsigned int i = 0; i < 3;i++){
		CI13[i] = ComputeIntersection<Simplex<1>, Simplex<3>>(triangle.getAbscissa(i) ,tetrahedron);
	}

};

void ComputeIntersection<Simplex<2>, Simplex<3>>::clear_all(){
	for(unsigned int i = 0; i < 6;i++){
		this->p_coordinates_tetrahedron[i] = NULL;
	}
	for(unsigned int i = 0; i < 3;i++){
		this->p_coordinates_triange[i] = NULL;
	}
};

void ComputeIntersection<Simplex<2>, Simplex<3>>::init(){
	for(unsigned int i = 0; i < 6;i++){
		CI12[i].setPC_abscissa(this->p_coordinates_tetrahedron[i]);
		CI12[i].setPC_triangle(this->p_coordinates_triange[0], 0);
		CI12[i].setPC_triangle(this->p_coordinates_triange[1], 1);
		CI12[i].setPC_triangle(this->p_coordinates_triange[2], 2);
	}
	for(unsigned int i = 0; i < 3;i++){
		CI13[i].setPC_abscissa(this->p_coordinates_triange[i]);
		CI13[i].setPC_tetrahedron(this->p_coordinates_tetrahedron[0],0);
		CI13[i].setPC_tetrahedron(this->p_coordinates_tetrahedron[1],1);
		CI13[i].setPC_tetrahedron(this->p_coordinates_tetrahedron[2],2);
		CI13[i].setPC_tetrahedron(this->p_coordinates_tetrahedron[3],3);
		CI13[i].setPC_tetrahedron(this->p_coordinates_tetrahedron[4],4);
		CI13[i].setPC_tetrahedron(this->p_coordinates_tetrahedron[5],5);
		CI13[i].init();
	}
};

void ComputeIntersection<Simplex<2>, Simplex<3>>::compute(){

};

Plucker* ComputeIntersection<Simplex<2>, Simplex<3>>::getPC_triangle(unsigned int index){
	return p_coordinates_triange[index];
};

Plucker* ComputeIntersection<Simplex<2>, Simplex<3>>::getPC_tetrahedron(unsigned int index){
	return p_coordinates_tetrahedron[index];
};

void ComputeIntersection<Simplex<2>, Simplex<3>>::setPC_triangle(Plucker* p_triangle_coordinate, unsigned int index){
	p_coordinates_triange[index] = p_triangle_coordinate;
};

void ComputeIntersection<Simplex<2>, Simplex<3>>::setPC_tetrahedron(Plucker* p_tetrahedron_coordinate, unsigned int index){
	p_coordinates_tetrahedron[index] = p_tetrahedron_coordinate;
};

} // END namespace_close
