/*
 *      Author: viktor
 */

#include "computeintersection.h"
#include "system/system.hh"
namespace computeintersection{

/****************************************************************
 * METODY PRO SIMPLEX 1 A SIMPLEX 2
 ****************************************************************/
ComputeIntersection<Simplex<1>, Simplex<2>>::ComputeIntersection(){
	abscissa = NULL;
	triangle = NULL;
};

ComputeIntersection<Simplex<1>, Simplex<2>>::ComputeIntersection(Simplex<1> &abs,Simplex<2> &triang){
	abscissa = &abs;
	triangle = &triang;

	plucker_coordinates_triangle.assign(3, new Plucker());
	plucker_coordinates_abscissa.assign(1, new Plucker);
	//plucker_coordinates_triangle.reserve(3);
	//plucker_coordinates_abscissa.reserve(1);
	//for(unsigned int i = 0; i < 3; i++){
			//plucker_coordinates_triangle[i] = new Plucker();
	//}
	//plucker_coordinates_abscissa[0] = new Plucker();
	this->clear_all();
};

void ComputeIntersection<Simplex<1>, Simplex<2>>::clear_all(){
	for(unsigned int i = 0; i < 3;i++){
		//this->p_coordinates_triangle[i] = NULL;
		plucker_products[i] = NULL;
	}
	//this->p_coordinates_abscissa[0] = NULL;
};

bool ComputeIntersection<Simplex<1>, Simplex<2>>::compute(double &theta, arma::vec3 &local_triangle){
	// Spočtení pluckerovych souradnic

	if(!plucker_coordinates_abscissa[0]->isComputed()){
		//cout << "pocitam abscissa" << endl;
		plucker_coordinates_abscissa[0]->compute((*abscissa)[0].getPointCoordinates(),
												 (*abscissa)[1].getPointCoordinates());
	}
	for(unsigned int i = 0; i < 3; i++){
		if(!plucker_coordinates_triangle[i]->isComputed()){
			//cout << "pocitam triangle" << endl;
			plucker_coordinates_triangle[i]->compute((*triangle)[i][0].getPointCoordinates(), (*triangle)[i][1].getPointCoordinates());
		}
	}

	// Vypočítání součinu dvou pluckerovych souřadnic a sdělení, jestli se jedná o průnik

	for(unsigned int i = 0; i < 3; i++){
		//cout << "[" << i << "]" ;
		if(plucker_products[i] == NULL){
		plucker_products[i] = new double((*plucker_coordinates_abscissa[0])*(*plucker_coordinates_triangle[i]));
		//cout << "PP nastaven = " << *plucker_products[i] << endl;
		}else{
		//cout << "PP pouzit = " << *plucker_products[i] << endl;
		}
	}

	if(((*plucker_products[0] > 0) && (*plucker_products[1] < 0) && (*plucker_products[2] > 0)) ||
	   ((*plucker_products[0] < 0) && (*plucker_products[1] > 0) && (*plucker_products[2] < 0))){
		double c = *plucker_products[0]; //c = -c;
		double d = *plucker_products[1]; d = -d;
		double e = *plucker_products[2]; //e = -e;

		// c = w0; d = w1; e = w2
		// lokální alfa = w2/soucet; lokální beta = w1/soucet; => lokální souřadnice na stěně
		local_triangle[0] = e/(c+d+e); // alfa
		local_triangle[1] = d/(c+d+e); // beta
		local_triangle[2] = c/(c+d+e); // gama

		// lokální souřadnice na přímce T
		// T = localAbscissa= (- A(i) + ( 1 - alfa - beta ) * V0(i) + alfa * V1(i) + beta * V2 (i)) / U(i)
		// i = max z U(i)
		arma::vec3 vec = (*abscissa)[1].getPointCoordinates() - (*abscissa)[0].getPointCoordinates();
		unsigned int i = 0;
		double max = vec[0];

		if(fabs(vec[1]) > fabs(max)){ max = vec[1]; i = 1;}
		if(fabs(vec[2]) > fabs(max)){ max = vec[2]; i = 2;}

		arma::vec3 global_triangle =
		local_triangle[0]*(*triangle)[0][0].getPointCoordinates() +
		local_triangle[1]*(*triangle)[0][1].getPointCoordinates() +
		local_triangle[2]*(*triangle)[1][1].getPointCoordinates();
		theta = (-(*abscissa)[0].getPointCoordinates()[i] + global_triangle[i])/max;
		//arma::vec3 local = ((-1)*abscissa[0].getPointCoordinates() + local_triangle)/vec;
		//arma::vec3 global_abscissa = local_abscissa * abscissa[1].getPointCoordinates() + (1 - local_abscissa) * abscissa[0].getPointCoordinates();
		xprintf(Msg, "vypocitany 3D (%f, %f, %f)\n", global_triangle[0], global_triangle[1], global_triangle[2]);
		return true;
	}else{
		return false;
	}
};

std::vector<Plucker *> &ComputeIntersection<Simplex<1>, Simplex<2>>::getPC_abscissa(){
	return plucker_coordinates_abscissa;
};
std::vector<Plucker *> &ComputeIntersection<Simplex<1>, Simplex<2>>::getPC_triangle(){
	return plucker_coordinates_triangle;
};
Plucker &ComputeIntersection<Simplex<1>, Simplex<2>>::getPC_abscissa(unsigned int index){
	return *plucker_coordinates_abscissa[index];
};
Plucker &ComputeIntersection<Simplex<1>, Simplex<2>>::getPC_triangle(unsigned int index){
	return *plucker_coordinates_triangle[index];
};

void ComputeIntersection<Simplex<1>, Simplex<2>>::setPC_abscissa(std::vector<Plucker *> &p_abscissa_coordinates){
	plucker_coordinates_abscissa = p_abscissa_coordinates;
};
void ComputeIntersection<Simplex<1>, Simplex<2>>::setPC_triangle(std::vector<Plucker *> &p_triangle_coordinates){
	plucker_coordinates_triangle= p_triangle_coordinates;
};
void ComputeIntersection<Simplex<1>, Simplex<2>>::setPC_abscissa(Plucker &p_abscissa_coordinate){
	plucker_coordinates_abscissa[0] = &p_abscissa_coordinate;
};
void ComputeIntersection<Simplex<1>, Simplex<2>>::setPC_triangle(Plucker &p_triangle_coordinate, unsigned int index){
	plucker_coordinates_triangle[index] = &p_triangle_coordinate;
};


void ComputeIntersection<Simplex<1>, Simplex<2>>::toStringPluckerCoordinates(){
	cout << "\tPluckerCoordinates Abscissa[0]";
		if(this->plucker_coordinates_abscissa[0] == NULL){
			cout << "NULL" << endl;
		}else{
			this->plucker_coordinates_abscissa[0]->toString();
		}
	for(unsigned int i = 0; i < 3;i++){
		cout << "\tPluckerCoordinates Triangle[" << i << "]";
		if(this->plucker_coordinates_triangle[i] == NULL){
			cout << "NULL" << endl;
		}else{
			this->plucker_coordinates_triangle[i]->toString();
		}
	}
};

void ComputeIntersection<Simplex<1>, Simplex<2>>::setPluckerProduct(double *number, unsigned int i){
	plucker_products[i] = number;
};

double* ComputeIntersection<Simplex<1>, Simplex<2>>::getPluckerProduct(unsigned int i){
	return plucker_products[i];
};


/****************************************************************
 * METODY PRO SIMPLEX 1 A SIMPLEX 3
 ****************************************************************/
ComputeIntersection<Simplex<1>, Simplex<3>>::ComputeIntersection(){
	abscissa = NULL;
	tetrahedron = NULL;
};

ComputeIntersection<Simplex<1>, Simplex<3>>::ComputeIntersection(Simplex<1> &abs,Simplex<3> &tetr){
	abscissa = &abs;
	tetrahedron = &tetr;

	plucker_coordinates_tetrahedron.assign(6, new Plucker());
	plucker_coordinates_abscissa.assign(1, new Plucker());

	for(unsigned int i = 0; i < 4;i++){
	 CI12[i] = ComputeIntersection<Simplex<1>, Simplex<2>>(*abscissa ,(*tetrahedron)[i]);
	}
};


void ComputeIntersection<Simplex<1>, Simplex<3>>::clear_all(){
	for(unsigned int i = 0; i < 6;i++){
		plucker_coordinates_tetrahedron[i]->clear();
	}
	plucker_coordinates_abscissa[0]->clear();

};

void ComputeIntersection<Simplex<1>, Simplex<3>>::init(){
	for(unsigned int i = 0; i < 4;i++){
		CI12[i].setPC_abscissa(*plucker_coordinates_abscissa[0]);
	 }
	CI12[0].setPC_triangle(*plucker_coordinates_tetrahedron[0], 0);
	CI12[0].setPC_triangle(*plucker_coordinates_tetrahedron[1], 1);
	CI12[0].setPC_triangle(*plucker_coordinates_tetrahedron[2], 2);
	CI12[1].setPC_triangle(*plucker_coordinates_tetrahedron[0], 0);
	CI12[1].setPC_triangle(*plucker_coordinates_tetrahedron[3], 1);
	CI12[1].setPC_triangle(*plucker_coordinates_tetrahedron[4], 2);
	CI12[2].setPC_triangle(*plucker_coordinates_tetrahedron[1], 0);
	CI12[2].setPC_triangle(*plucker_coordinates_tetrahedron[3], 1);
	CI12[2].setPC_triangle(*plucker_coordinates_tetrahedron[5], 2);
	CI12[3].setPC_triangle(*plucker_coordinates_tetrahedron[2], 0);
	CI12[3].setPC_triangle(*plucker_coordinates_tetrahedron[4], 1);
	CI12[3].setPC_triangle(*plucker_coordinates_tetrahedron[5], 2);

};

void ComputeIntersection<Simplex<1>, Simplex<3>>::compute(IntersectionLocal &lokalni_mnohouhelnik, unsigned int index_edge = 0){

	//xprintf(Msg, "Pruchod metodou compute 13\n");

	//std::vector<double> docasna;docasna.assign(3,0);
	arma::vec3 docasna;docasna[0] = 0;docasna[1] = 0;docasna[2] = 0;
	double c = 0;
	unsigned int pocet_pruniku = 0;

	if(CI12[0].compute(c, docasna)){

		IntersectionPoint<1,2> IP(c, docasna,index_edge, 0);
		lokalni_mnohouhelnik.add_local_point(IntersectionLocal::convertTo23(IP));
		pocet_pruniku++;
	}

	CI12[1].setPluckerProduct(CI12[0].getPluckerProduct(0),0);

	if(CI12[1].compute(c, docasna)){
		docasna[2] = 0;
		IntersectionPoint<1,2> IP(c, docasna,index_edge,1);
		lokalni_mnohouhelnik.add_local_point(IntersectionLocal::convertTo23(IP));
				arma::vec3 bod = docasna[0]*(tetrahedron->getAbscissa(0))[0].getPointCoordinates() +
						docasna[1]*(tetrahedron->getAbscissa(0))[1].getPointCoordinates() +
						docasna[2]*(tetrahedron->getAbscissa(1))[1].getPointCoordinates() +
						(1 - docasna[0] - docasna[1] - docasna[2])*(tetrahedron->getAbscissa(3))[1].getPointCoordinates();
				xprintf(Msg, "Globální průnik: (%f,%f,%f) \n", bod[0], bod[1], bod[2]);
				pocet_pruniku++;
	}

	if(pocet_pruniku >= 2){
		return;
	}

	CI12[2].setPluckerProduct(CI12[0].getPluckerProduct(1),0);
	CI12[2].setPluckerProduct(CI12[1].getPluckerProduct(1),1);
	if(CI12[2].compute(c, docasna)){
		docasna[2] = docasna[1];
		docasna[1] = 0;
		IntersectionPoint<1,2> IP(c, docasna,index_edge,2);
		lokalni_mnohouhelnik.add_local_point(IntersectionLocal::convertTo23(IP));
				arma::vec3 bod = docasna[0]*(tetrahedron->getAbscissa(0))[0].getPointCoordinates() +
						docasna[1]*(tetrahedron->getAbscissa(0))[1].getPointCoordinates() +
						docasna[2]*(tetrahedron->getAbscissa(1))[1].getPointCoordinates() +
						(1 - docasna[0] - docasna[1] - docasna[2])*(tetrahedron->getAbscissa(3))[1].getPointCoordinates();
				xprintf(Msg, "Globální průnik: (%f,%f,%f) \n", bod[0], bod[1], bod[2]);
				pocet_pruniku++;
	}

	if(pocet_pruniku >= 2){
			return;
	}

	CI12[3].setPluckerProduct(CI12[0].getPluckerProduct(2),0);
	CI12[3].setPluckerProduct(CI12[1].getPluckerProduct(2),1);
	CI12[3].setPluckerProduct(CI12[2].getPluckerProduct(2),2);
	if(CI12[3].compute(c, docasna)){
		docasna[2] = docasna[1];
		docasna[1] = docasna[0];
		docasna[0] = 0;
		IntersectionPoint<1,2> IP(c, docasna,index_edge,3);
		lokalni_mnohouhelnik.add_local_point(IntersectionLocal::convertTo23(IP));
				arma::vec3 bod = docasna[0]*(tetrahedron->getAbscissa(0))[0].getPointCoordinates() +
						docasna[1]*(tetrahedron->getAbscissa(0))[1].getPointCoordinates() +
						docasna[2]*(tetrahedron->getAbscissa(1))[1].getPointCoordinates() +
						(1 - docasna[0] - docasna[1] - docasna[2])*(tetrahedron->getAbscissa(3))[1].getPointCoordinates();
				xprintf(Msg, "Globální průnik: (%f,%f,%f) \n", bod[0], bod[1], bod[2]);
				pocet_pruniku++;
	}
};

std::vector<Plucker *> &ComputeIntersection<Simplex<1>, Simplex<3>>::getPC_abscissa(){
	return plucker_coordinates_abscissa;
};
std::vector<Plucker *> &ComputeIntersection<Simplex<1>, Simplex<3>>::getPC_tetrahedron(){
	return plucker_coordinates_tetrahedron;
};
Plucker &ComputeIntersection<Simplex<1>, Simplex<3>>::getPC_abscissa(unsigned int index){
	return *plucker_coordinates_abscissa[index];
};
Plucker &ComputeIntersection<Simplex<1>, Simplex<3>>::getPC_tetrahedron(unsigned int index){
	return *plucker_coordinates_tetrahedron[index];
};

void ComputeIntersection<Simplex<1>, Simplex<3>>::setPC_abscissa(std::vector<Plucker *> &p_abscissa_coordinates){
	plucker_coordinates_abscissa = p_abscissa_coordinates;
};
void ComputeIntersection<Simplex<1>, Simplex<3>>::setPC_tetrahedron(std::vector<Plucker *> &p_tetrahedron_coordinates){
	plucker_coordinates_tetrahedron = p_tetrahedron_coordinates;
};
void ComputeIntersection<Simplex<1>, Simplex<3>>::setPC_abscissa(Plucker &p_abscissa_coordinate){
	plucker_coordinates_abscissa[0] = &p_abscissa_coordinate;
};
void ComputeIntersection<Simplex<1>, Simplex<3>>::setPC_tetrahedron(Plucker &p_tetrahedron_coordinate, unsigned int index){
	plucker_coordinates_tetrahedron[index] = &p_tetrahedron_coordinate;
};



void ComputeIntersection<Simplex<1>, Simplex<3>>::toStringPluckerCoordinates(){
		cout << "\tPluckerCoordinates Abscissa[0]";
		if(plucker_coordinates_abscissa[0] == NULL){
			cout << "NULL" << endl;
		}else{
			this->plucker_coordinates_abscissa[0]->toString();
		}

	for(unsigned int i = 0; i < 6;i++){
		cout << "\tPluckerCoordinates Tetrahedron[" << i << "]";
		if(plucker_coordinates_tetrahedron[i] == NULL){
			cout << "NULL" << endl;
		}else{
			this->plucker_coordinates_tetrahedron[i]->toString();
		}
	}
};

void ComputeIntersection<Simplex<1>, Simplex<3>>::toStringPluckerCoordinatesTree(){
	cout << "ComputeIntersection<Simplex<1>, <Simplex<3>> Plucker Coordinates Tree:" << endl;
		this->toStringPluckerCoordinates();
		for(unsigned int i = 0; i < 4;i++){
			cout << "ComputeIntersection<Simplex<1>, Simplex<2>>["<< i <<"] Plucker Coordinates:" << endl;
			CI12[i].toStringPluckerCoordinates();
		}
};

void ComputeIntersection<Simplex<1>, Simplex<3>>::setPluckerProduct(double* number, unsigned int index_CI, unsigned index_edge){
	CI12[index_CI].setPluckerProduct(number, index_edge);
};

/****************************************************************
 * METODY PRO SIMPLEX 2 A SIMPLEX 3
 ****************************************************************/
ComputeIntersection<Simplex<2>, Simplex<3>>::ComputeIntersection(){
	this->triange = NULL;
	this->tetrahedron = NULL;
};


ComputeIntersection<Simplex<2>, Simplex<3>>::ComputeIntersection(Simplex<2> &triangle, Simplex<3> &tetr){
	this->triange = &triangle;
	this->tetrahedron = &tetr;

	plucker_coordinates_triangle.reserve(3);
	plucker_coordinates_tetrahedron.reserve(6);

	for(unsigned int i = 0; i < 6;i++){
		plucker_coordinates_tetrahedron[i] = new Plucker();
		CI12[i] = ComputeIntersection<Simplex<1>, Simplex<2>>(tetrahedron->getAbscissa(i) ,*triange);
	}
	for(unsigned int i = 0; i < 3;i++){
		plucker_coordinates_triangle[i] = new Plucker();
		CI13[i] = ComputeIntersection<Simplex<1>, Simplex<3>>(triange->getAbscissa(i) , *tetrahedron);
	}
};

void ComputeIntersection<Simplex<2>, Simplex<3>>::clear_all(){
	for(unsigned int i = 0; i < 3;i++){
		plucker_coordinates_triangle[i]->clear();
		plucker_coordinates_tetrahedron[2*i]->clear();
		plucker_coordinates_tetrahedron[2*i+1]->clear();
	}
};

void ComputeIntersection<Simplex<2>, Simplex<3>>::init(){


for(unsigned int i = 0; i < 6;i++){
	CI12[i].setPC_abscissa(*plucker_coordinates_tetrahedron[i]);
	//CI12[i].setPC_triangle(plucker_coordinates_triangle);
	CI12[i].setPC_triangle(*plucker_coordinates_triangle[0], 0);
	CI12[i].setPC_triangle(*plucker_coordinates_triangle[1], 1);
	CI12[i].setPC_triangle(*plucker_coordinates_triangle[2], 2);
	}
	for(unsigned int i = 0; i < 3;i++){
		CI13[i].setPC_abscissa(*plucker_coordinates_triangle[i]);
		CI13[i].setPC_tetrahedron(*plucker_coordinates_tetrahedron[0],0);
		CI13[i].setPC_tetrahedron(*plucker_coordinates_tetrahedron[1],1);
		CI13[i].setPC_tetrahedron(*plucker_coordinates_tetrahedron[2],2);
		CI13[i].setPC_tetrahedron(*plucker_coordinates_tetrahedron[3],3);
		CI13[i].setPC_tetrahedron(*plucker_coordinates_tetrahedron[4],4);
		CI13[i].setPC_tetrahedron(*plucker_coordinates_tetrahedron[5],5);
		//CI13[i].setPC_tetrahedron(plucker_coordinates_tetrahedron);
		CI13[i].init();
		/*CI13[i].setPC_abscissa(this->p_coordinates_triange[i]);

		CI13[i].init();*/
	}
};

void ComputeIntersection<Simplex<2>, Simplex<3>>::compute(IntersectionLocal &lokalni_mnohouhlenik){

	// čekat dokud se nenalezne první průsečík -> pak pokračovat jiným algoritmem
	// == tady bude optimalizovaný algoritmus




	// hrubý algoritmus
	arma::vec3 docasna;
	double c = 0;

	for(unsigned int i = 0; i < 6;i++){
		CI12[i].compute(c, docasna);
	}

	// Optimalizace: znovu použití již vypočítaných součinů
	for(unsigned int i = 0; i < 3;i++){
		CI13[i].setPluckerProduct(CI12[0].getPluckerProduct(i),0,0);
		CI13[i].setPluckerProduct(CI12[1].getPluckerProduct(i),0,1);
		CI13[i].setPluckerProduct(CI12[2].getPluckerProduct(i),0,2);
		CI13[i].setPluckerProduct(CI12[3].getPluckerProduct(i),1,1);
		CI13[i].setPluckerProduct(CI12[4].getPluckerProduct(i),1,2);
		CI13[i].setPluckerProduct(CI12[5].getPluckerProduct(i),2,2);
	}


	for(unsigned int i = 0; i < 3;i++){
		CI13[i].compute(lokalni_mnohouhlenik);
	}
};



	std::vector<Plucker *> &ComputeIntersection<Simplex<2>, Simplex<3>>::getPC_triangle(){
		return plucker_coordinates_triangle;
	};
	std::vector<Plucker *> &ComputeIntersection<Simplex<2>, Simplex<3>>::getPC_tetrahedron(){
		return plucker_coordinates_tetrahedron;
	};
	Plucker &ComputeIntersection<Simplex<2>, Simplex<3>>::getPC_triangle(unsigned int index){
		return *plucker_coordinates_triangle[index];
	};
	Plucker &ComputeIntersection<Simplex<2>, Simplex<3>>::getPC_tetrahedron(unsigned int index){
		return *plucker_coordinates_tetrahedron[index];
	};

	void ComputeIntersection<Simplex<2>, Simplex<3>>::setPC_triangle(std::vector<Plucker *> &p_triangle_coordinates){
		plucker_coordinates_triangle = p_triangle_coordinates;
	};
	void ComputeIntersection<Simplex<2>, Simplex<3>>::setPC_tetrahedron(std::vector<Plucker *> &p_tetrahedron_coordinates){
		plucker_coordinates_tetrahedron = p_tetrahedron_coordinates;
	};
	void ComputeIntersection<Simplex<2>, Simplex<3>>::setPC_triangle(Plucker &p_triangle_coordinate, unsigned int index){
		plucker_coordinates_triangle[index] = &p_triangle_coordinate;
	};
	void ComputeIntersection<Simplex<2>, Simplex<3>>::setPC_tetrahedron(Plucker &p_tetrahedron_coordinate, unsigned int index){
		plucker_coordinates_tetrahedron[index] = &p_tetrahedron_coordinate;
	};


void ComputeIntersection<Simplex<2>, Simplex<3>>::toStringPluckerCoordinates(){
	for(unsigned int i = 0; i < 3;i++){
		cout << "\tPluckerCoordinates Triangle[" << i << "]";
		if(plucker_coordinates_triangle[i] == NULL){
			cout << "NULL" << endl;
		}else{
			plucker_coordinates_triangle[i]->toString();
		}
	}
	for(unsigned int i = 0; i < 6;i++){
		cout << "\tPluckerCoordinates Tetrahedron[" << i << "]";
		if(plucker_coordinates_tetrahedron[i] == NULL){
			cout << "NULL" << endl;
		}else{
			plucker_coordinates_tetrahedron[i]->toString();
		}
	}
};

void ComputeIntersection<Simplex<2>, Simplex<3>>::toStringPluckerCoordinatesTree(){
	cout << "ComputeIntersection<Simplex<2>, <Simplex<3>> Plucker Coordinates Tree:" << endl;
	this->toStringPluckerCoordinates();
	for(unsigned int i = 0; i < 6;i++){
		cout << "ComputeIntersection<Simplex<1>, Simplex<2>>["<< i <<"] Plucker Coordinates:" << endl;
		CI12[i].toStringPluckerCoordinates();
	}
	for(unsigned int i = 0; i < 3;i++){
		CI13[i].toStringPluckerCoordinatesTree();
	}
};

} // END namespace_close
