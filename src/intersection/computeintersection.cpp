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
	computed = false;
	abscissa = NULL;
	triangle = NULL;
	plucker_coordinates_triangle.reserve(3);
	plucker_coordinates_abscissa.reserve(1);
};

ComputeIntersection<Simplex<1>, Simplex<2>>::ComputeIntersection(Simplex<1> &abs,Simplex<2> &triang){
	computed = false;
	abscissa = &abs;
	triangle = &triang;

	plucker_coordinates_triangle.assign(3, new Plucker());
	plucker_coordinates_abscissa.assign(1, new Plucker);

	this->clear_all();
};

void ComputeIntersection<Simplex<1>, Simplex<2>>::clear_all(){
	for(unsigned int i = 0; i < 3;i++){
		//this->p_coordinates_triangle[i] = NULL;
		plucker_products[i] = NULL;
	}
	//this->p_coordinates_abscissa[0] = NULL;
};

void ComputeIntersection<Simplex<1>, Simplex<2>>::initPluckerToCompute(){
	// Spočtení pluckerovych souradnic

	if(!plucker_coordinates_abscissa[0]->is_computed()){
		plucker_coordinates_abscissa[0]->compute((*abscissa)[0].getPointCoordinates(),
												 (*abscissa)[1].getPointCoordinates());
	}
	for(unsigned int i = 0; i < 3; i++){
		if(!plucker_coordinates_triangle[i]->is_computed()){
			plucker_coordinates_triangle[i]->compute((*triangle)[i][0].getPointCoordinates(), (*triangle)[i][1].getPointCoordinates());
		}
	}

	for(unsigned int i = 0; i < 3; i++){
		if(plucker_products[i] == NULL){
			plucker_products[i] = new double((*plucker_coordinates_abscissa[0])*(*plucker_coordinates_triangle[i]));
		}
	}

};

void ComputeIntersection<Simplex<1>, Simplex<2>>::set_data(Simplex<1> *abs, Simplex<2> *triang){
	computed = false;
	abscissa = abs;
	triangle = triang;
	this->clear_all();
};

bool ComputeIntersection<Simplex<1>, Simplex<2>>::isComputed(){
	return computed;
};

void ComputeIntersection<Simplex<1>, Simplex<2>>::setComputed(){
	computed = true;
};

bool ComputeIntersection<Simplex<1>, Simplex<2>>::compute(std::vector<IntersectionPoint<1,2>> &IP12s, bool compute_zeros_plucker_products){

	initPluckerToCompute();

	computed = true;

	if(((*plucker_products[0] > 0) && (*plucker_products[1] < 0) && (*plucker_products[2] > 0)) ||
	   ((*plucker_products[0] < 0) && (*plucker_products[1] > 0) && (*plucker_products[2] < 0))){
		double c = *plucker_products[0]; //c = -c;
		double d = *plucker_products[1]; d = -d;
		double e = *plucker_products[2]; //e = -e;
		//xprintf(Msg,"Prunik.%f %f %f\n",c,d,e);
		// c = w0; d = w1; e = w2
		// lokální alfa = w2/soucet; lokální beta = w1/soucet; => lokální souřadnice na stěně
		arma::vec::fixed<3> local_triangle;
		arma::vec::fixed<2> theta;
		local_triangle[0] = e/(c+d+e); // alfa
		local_triangle[1] = d/(c+d+e); // beta
		local_triangle[2] = c/(c+d+e); // gama

		// lokální souřadnice na přímce T
		// T = localAbscissa= (- A(i) + ( 1 - alfa - beta ) * V0(i) + alfa * V1(i) + beta * V2 (i)) / U(i)
		// i = max z U(i)
		arma::vec3 vec = (*abscissa)[1].getPointCoordinates() - (*abscissa)[0].getPointCoordinates();
		unsigned int i = 0;
		double max = vec[0];

		if(fabs((double)vec[1]) > fabs(max)){ max = vec[1]; i = 1;}
		if(fabs((double)vec[2]) > fabs(max)){ max = vec[2]; i = 2;}

		arma::vec3 global_triangle =
		local_triangle[0]*(*triangle)[0][0].getPointCoordinates() +
		local_triangle[1]*(*triangle)[0][1].getPointCoordinates() +
		local_triangle[2]*(*triangle)[1][1].getPointCoordinates();
		theta[1] = (-(*abscissa)[0].getPointCoordinates()[i] + global_triangle[i])/max;
		theta[0] = 1 - theta[1];


		IntersectionPoint<1,2> IP(theta,local_triangle,-1,-1,(*plucker_products[0] > 0 ? 1 : 0));
		//IP12s[0] = IP;
		IP12s.push_back(IP);
		return true;
	}else if(compute_zeros_plucker_products && (*plucker_products[0] == 0 || *plucker_products[1] == 0 || *plucker_products[2] == 0)){

		unsigned int num_zeros = 0;
		for(unsigned int i = 0; i < 3;i++){
			if(*plucker_products[i] == 0){
				num_zeros++;
			}
		}

		// Pokud je 1 nula = jeden patologický průsečík

		// Pokud jsou 2 nuly => průsečík ve vrcholu

		// Pokud jsou 3 nuly - všechny vypočítat

		for(unsigned int i = 0; i < 3;i++){
			if(*plucker_products[i] == 0){
				arma::vec3 A = (*abscissa)[0].getPointCoordinates();
				//arma::vec3 B("9 5 0");
				arma::vec3 U = plucker_coordinates_abscissa[0]->get_u_vector();
				arma::vec3 C = (*triangle)[i][i%2].getPointCoordinates();
				//arma::vec3 D("4 4 0");
				arma::vec3 V = plucker_coordinates_triangle[i]->get_u_vector();
				arma::vec3 K = C - A;
				arma::vec3 Det = -arma::cross(U,V);
				unsigned int max_index = 0;
				double maximum = Det[0];
				if(fabs((double)Det[1]) > fabs(maximum)){
					maximum = Det[1];
					max_index = 1;
				}
				if(fabs((double)Det[2]) > fabs(maximum)){
					maximum = Det[2];
					max_index = 2;
				}
				if(maximum == 0){
					return false;
				}

				double DetX = K[(max_index+2)%3]*V[(max_index+1)%3]
							-K[(max_index+1)%3]*V[(max_index+2)%3];

			    double DetY = K[(max_index+2)%3]*U[(max_index+1)%3]
							-K[(max_index+1)%3]*U[(max_index+2)%3];

				double s;
				double t;

				s = DetX/Det[max_index];
				t = DetY/Det[max_index];

				if(i == 1){
					t = -t;
				}

				//xprintf(Msg, "s,t : %f,%f \n",s,t);

				if(t > 1 || t < 0){
					//return false;
				}else{

					arma::vec::fixed<2> local_abscissa;
					local_abscissa[0] = 1-s;
					local_abscissa[1] = s;

					arma::vec::fixed<3> l_triangle;
					/*
					 * 0 = 1- t, t, 0
					 * 1 = 1 -t , 0, t // 1 = t, 0 , 1 - t
					 * 2 = 0, 1 - t, t
					 *
					 * t = 0 => 1; 1 => 0, 2 => 2
					 *
					 * */
					l_triangle[(3-i)%3] = 1 - t;
					l_triangle[(4-i)%3] = t;
					l_triangle[2-i] = 0;

					IntersectionPoint<1,2> IP(local_abscissa,l_triangle,-1,i,1,false,true);
					IP12s.push_back(IP);
					return true;
				}




			}
		}
		return false;
	}else{
		return false;
	}
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
	plucker_coordinates_tetrahedron.reserve(6);
	plucker_coordinates_abscissa.reserve(1);
};

ComputeIntersection<Simplex<1>, Simplex<3>>::ComputeIntersection(Simplex<1> &abs,Simplex<3> &tetr){
	abscissa = &abs;
	tetrahedron = &tetr;

	plucker_coordinates_tetrahedron.assign(6, new Plucker());
	plucker_coordinates_abscissa.assign(1, new Plucker());

	for(unsigned int i = 0; i < 4;i++){
		CI12[i].set_data(abscissa, &(*tetrahedron)[i]);
	 //CI12[i] = ComputeIntersection<Simplex<1>, Simplex<2>>(*abscissa ,(*tetrahedron)[i]);
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
		CI12[i].set_pc_abscissa(plucker_coordinates_abscissa[0]);
	}

	for(unsigned int j = 0; j < 4;j++){
		for(unsigned int i = 0; i < 3;i++){
			CI12[j].set_pc_triangle(plucker_coordinates_tetrahedron[RefSimplex<3>::side_lines[j][i]], i);
		}
	}

};

void ComputeIntersection<Simplex<1>, Simplex<3>>::set_data(Simplex<1> *abs, Simplex<3> *tetr){
	abscissa = abs;
	tetrahedron = tetr;
	for(unsigned int i = 0; i < 4;i++){
		CI12[i].set_data(abscissa, &(*tetrahedron)[i]);
	}
};

int ComputeIntersection<Simplex<1>, Simplex<3>>::compute(std::vector<IntersectionPoint<1,3>> &IP13s){

	std::vector<IntersectionPoint<1,2>> IP12s;
	unsigned int pocet_pruniku = 0;

	for(unsigned int i = 0;i < 4 && pocet_pruniku < 2;i++){

		for(unsigned int j = 0; j < i;j++){
			CI12[i].setPluckerProduct(CI12[j].getPluckerProduct(i-1),j);
		}
		if(!CI12[i].isComputed() && CI12[i].compute(IP12s, true)){
			//xprintf(Msg,"Prunik13\n");

			if(IP12s[IP12s.size() - 1].isPatological()){
				//xprintf(Msg, "\tPatologicky\n");
				// Nastavování stěn, které se už nemusí počítat
				if(i == 0){
					if(IP12s[IP12s.size() - 1].getLocalCoords2()[0] == 1){
						CI12[1].setComputed();
						CI12[2].setComputed();
					}else if(IP12s[IP12s.size() - 1].getLocalCoords2()[1] == 1){
						CI12[1].setComputed();
						CI12[3].setComputed();
					}else if(IP12s[IP12s.size() - 1].getLocalCoords2()[2] == 1){
						CI12[2].setComputed();
						CI12[3].setComputed();
					}else{
						CI12[IP12s[IP12s.size() - 1].getSide2() + 1].setComputed();
					}
				}else if(i == 1 && IP12s[IP12s.size() - 1].getLocalCoords2()[2] == 1){
					CI12[2].setComputed();
					CI12[3].setComputed();
				}else{
					CI12[IP12s[IP12s.size() - 1].getSide2() + 1].setComputed();
				}
				IP12s[IP12s.size() - 1].setSide2(i);
			}else{
				IP12s[IP12s.size() - 1].setSide2(i);
			}
			pocet_pruniku++;

			//if((IP.getLocalCoords1())[0] <= 1 && (IP.getLocalCoords1())[0] >= 0){
				//IP.print();
				//xprintf(Msg, "Interpolace dimenzi");
				IntersectionPoint<1,3> IP13 = IntersectionLocal::interpolateDimension<1,3>(IP12s[IP12s.size() - 1]);
				//IP13.print();
				IP13s.push_back(IP13);

			//}
		}
	}

	// Kontrola vytvořených průniků => zda-li není potřeba interpolovat + zda-li se o prunik vubec nejedna:
	if(pocet_pruniku > 1){
		double first_theta = IP13s[IP13s.size()-2].getLocalCoords1()[1];
		double second_theta = IP13s[IP13s.size()-1].getLocalCoords1()[1];

		  // Nejedná se o průnik - celá usečka leží mimo čtyřstěn
		if(((first_theta > 1) && (second_theta > 1)) ||
		   ((first_theta < 0) && (second_theta < 0))){

			pocet_pruniku = 0;
			IP13s.pop_back();
			IP13s.pop_back();
		}else{

			// Jedná se o průnik
			// První souřadnice leží uvnitř čtyřstěnu
			if(first_theta > 1 || first_theta < 0){
				double theta = first_theta > 1 ? 1 : 0;
				arma::vec::fixed<4> interpolovane = RefSimplex<3>::line_barycentric_interpolation(IP13s[IP13s.size()-2].getLocalCoords2(), IP13s[IP13s.size()-1].getLocalCoords2(), first_theta, second_theta,theta);
				arma::vec::fixed<2> inter; inter[0] = 1 - theta; inter[1] = theta;
				IP13s[IP13s.size()-2].setLocalCoords2(interpolovane);
				IP13s[IP13s.size()-2].setLocalCoords1(inter);

					IP13s[IP13s.size()-2].setIsVertex(true);

				first_theta = theta;
			}
			// Druhá souřadnice leží uvnitř čtyřstěnu
			if(second_theta > 1 || second_theta < 0){
				double theta2 = second_theta > 1 ? 1 : 0;
				arma::vec::fixed<2> inter2; inter2[0] = 1 - theta2; inter2[1] = theta2;
				arma::vec::fixed<4> interpolovane2 = RefSimplex<3>::line_barycentric_interpolation(IP13s[IP13s.size()-2].getLocalCoords2(), IP13s[IP13s.size()-1].getLocalCoords2(), first_theta, second_theta,theta2);
				IP13s[IP13s.size()-1].setLocalCoords2(interpolovane2);
				IP13s[IP13s.size()-1].setLocalCoords1(inter2);

					IP13s[IP13s.size()-1].setIsVertex(true);

			}
		}
	}
	return pocet_pruniku;

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

double* ComputeIntersection<Simplex<1>, Simplex<3>>::getPluckerProduct(unsigned int index_CI, unsigned index_edge){
	return CI12[index_CI].getPluckerProduct(index_edge);
};

/****************************************************************
 * METODY PRO SIMPLEX 2 A SIMPLEX 3
 ****************************************************************/
ComputeIntersection<Simplex<2>, Simplex<3>>::ComputeIntersection(){
	this->triange = NULL;
	this->tetrahedron = NULL;

	plucker_coordinates_triangle.reserve(3);
	plucker_coordinates_tetrahedron.reserve(6);
};


ComputeIntersection<Simplex<2>, Simplex<3>>::ComputeIntersection(Simplex<2> &triangle, Simplex<3> &tetr){
	this->triange = &triangle;
	this->tetrahedron = &tetr;

	plucker_coordinates_triangle.reserve(3);
	plucker_coordinates_tetrahedron.reserve(6);

	for(unsigned int i = 0; i < 6;i++){
		plucker_coordinates_tetrahedron[i] = new Plucker();
		CI12[i].set_data(&tetrahedron->getAbscissa(i), triange);
	}
	for(unsigned int i = 0; i < 3;i++){
		plucker_coordinates_triangle[i] = new Plucker();
		CI13[i].set_data(&triange->getAbscissa(i) , tetrahedron);
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

		CI12[i].set_pc_abscissa(plucker_coordinates_tetrahedron[i]);
		for(unsigned int j = 0; j < 3;j++){
			CI12[i].set_pc_triangle(plucker_coordinates_triangle[j],j);
		}
	}
	for(unsigned int i = 0; i < 3;i++){
		CI13[i].set_pc_abscissa(plucker_coordinates_triangle[i]);
		for(unsigned int j = 0; j < 6;j++){
			CI13[i].set_pc_tetrahedron(plucker_coordinates_tetrahedron[j],j);
		}
		CI13[i].init();
	}
};

void ComputeIntersection<Simplex<2>, Simplex<3>>::compute(IntersectionLocal &lokalni_mnohouhelnik){

	std::vector<IntersectionPoint<1,2>> IP12s;
	std::vector<IntersectionPoint<1,3>> IP13s;
	int pocet_pruniku = 0;
	int pocet_13_pruniku;

	//cout << "ComputeIntersection<Simplex<2>, Simplex<3>>::compute - edges triangle vs tetrahedron" << endl;
		for(unsigned int i = 0; i < 3;i++){
			pocet_13_pruniku = CI13[(3-i)%3].compute(IP13s);
			//(triange->getAbscissa(i)).toString();
			// Vždy by měl být počet průniku 2 nebo 0
			if(pocet_13_pruniku == 1){
				IP13s[IP13s.size() - 1].setSide1((3-i)%3);
				IntersectionPoint<3,1> IP31 = IntersectionLocal::flipDimension<3,1>(IP13s[IP13s.size() - 1]);
				IntersectionPoint<3,2> IP32 = IntersectionLocal::interpolateDimension<3,2>(IP31);
				IntersectionPoint<2,3> IP23 = IntersectionLocal::flipDimension<2,3>(IP32);
				lokalni_mnohouhelnik.addIP(IP23);

			}else if(pocet_13_pruniku == 2){

				IP13s[IP13s.size() - 2].setSide1((3-i)%3);
				IntersectionPoint<3,1> IP31 = IntersectionLocal::flipDimension<3,1>(IP13s[IP13s.size() - 2]);
				IntersectionPoint<3,2> IP32 = IntersectionLocal::interpolateDimension<3,2>(IP31);
				IntersectionPoint<2,3> IP23 = IntersectionLocal::flipDimension<2,3>(IP32);
				lokalni_mnohouhelnik.addIP(IP23);

				IP13s[IP13s.size() - 1].setSide1((3-i)%3);
				IP31 = IntersectionLocal::flipDimension<3,1>(IP13s[IP13s.size() - 1]);
				IP32 = IntersectionLocal::interpolateDimension<3,2>(IP31);
				IP23 = IntersectionLocal::flipDimension<2,3>(IP32);
				lokalni_mnohouhelnik.addIP(IP23);

			}


		}
	// Optimalizace: znovu použití již vypočítaných součinů


	for(unsigned int i = 0; i < 3;i++){
		CI12[0].setPluckerProduct(CI13[i].getPluckerProduct(0,0),i);
		CI12[1].setPluckerProduct(CI13[i].getPluckerProduct(0,1),i);
		CI12[2].setPluckerProduct(CI13[i].getPluckerProduct(0,2),i);
		CI12[3].setPluckerProduct(CI13[i].getPluckerProduct(1,1),i);
		CI12[4].setPluckerProduct(CI13[i].getPluckerProduct(1,2),i);
		CI12[5].setPluckerProduct(CI13[i].getPluckerProduct(2,2),i);
	}

	//cout << "ComputeIntersection<Simplex<2>, Simplex<3>>::compute - edges tetrahedron vs triangle" << endl;
		for(unsigned int i = 0; i < 6;i++){
			if(CI12[i].compute(IP12s, false)){
				//xprintf(Msg,"Prunik21 %d\n", IP12s.size());
				pocet_pruniku++;
				IP12s[IP12s.size() - 1].setSide1(i);
				if((IP12s[IP12s.size() - 1].getLocalCoords1())[0] <= 1 && (IP12s[IP12s.size() - 1].getLocalCoords1())[0] >= 0){
					if(IP12s[IP12s.size() - 1].getLocalCoords1()[0] == 1 || IP12s[IP12s.size() - 1].getLocalCoords1()[0] == 0){
						IP12s[IP12s.size() - 1].setIsPatological(true);
					}
					//IP.print();
									IntersectionPoint<2,1> IP21 = IntersectionLocal::flipDimension<2,1>(IP12s[IP12s.size() - 1]);
									IntersectionPoint<2,3> IP23 = IntersectionLocal::interpolateDimension<2,3>(IP21);
									//IntersectionPoint<2,3> IP23 = IntersectionLocal::interpolateDimension<2,2>(IP22);
									//IP23.print();
									lokalni_mnohouhelnik.addIP(IP23);


				}

				//IP12s.pop_back();
			}
		}


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
