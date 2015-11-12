/*
 *      Author: viktor
 */

#include "computeintersection.h"
#include "mesh/ref_element.hh"
#include "system/system.hh"

namespace computeintersection{

/****************************************************************
 * METODY PRO SIMPLEX 1 A SIMPLEX 2
 ****************************************************************/
const double ComputeIntersection<Simplex<1>, Simplex<2>>::epsilon = 0.0000001;//128*2048*numeric_limits<double>::epsilon();

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

	plucker_coordinates_abscissa.reserve(1);
	plucker_coordinates_triangle.reserve(3);

	plucker_coordinates_abscissa[0] = new Plucker();

	for(unsigned int i = 0; i < 3;i++){
		plucker_coordinates_triangle[i] = new Plucker();
	}

	//plucker_coordinates_triangle.assign(3, new Plucker());
	//plucker_coordinates_abscissa.assign(1, new Plucker);

	clear_all();
};

void ComputeIntersection<Simplex<1>, Simplex<2>>::clear_all(){
	for(unsigned int i = 0; i < 3;i++){
		//p_coordinates_triangle[i] = NULL;
		plucker_products[i] = NULL;
	}
	//p_coordinates_abscissa[0] = NULL;
};

void ComputeIntersection<Simplex<1>, Simplex<2>>::init_plucker_to_compute(){

    // if not already computed, compute plucker coordinates of abscissa
	if(!plucker_coordinates_abscissa[0]->is_computed()){
		plucker_coordinates_abscissa[0]->compute((*abscissa)[0].getPointCoordinates(),
												 (*abscissa)[1].getPointCoordinates());
	}
	// if not already computed, compute plucker coordinates of triangle sides
	for(unsigned int side = 0; side < RefElement<2>::n_sides; side++){
		if(!plucker_coordinates_triangle[side]->is_computed()){
			plucker_coordinates_triangle[side]->compute((*triangle)[side][0].getPointCoordinates(), (*triangle)[side][1].getPointCoordinates());
		}
	}

// 	DBGMSG("Abscissa:\n");
//     (*abscissa)[0].getPointCoordinates().print();
//     (*abscissa)[1].getPointCoordinates().print();
    
	// compute Plucker products (abscissa X triangle side)
	for(unsigned int side = 0; side < RefElement<2>::n_sides; side++){
		if(plucker_products[side] == NULL){
//             (*plucker_coordinates_abscissa[0]).toString();
//             (*plucker_coordinates_triangle[side]).toString();
			plucker_products[side] = new double((*plucker_coordinates_abscissa[0])*(*plucker_coordinates_triangle[side]));
            
//             DBGMSG("triangle side:\n");
//             (*triangle)[side][0].getPointCoordinates().print();
//             (*triangle)[side][1].getPointCoordinates().print();
		}
// 		DBGMSG("Plucker product = %f\n", *(plucker_products[side]));
	}

};

void ComputeIntersection<Simplex<1>, Simplex<2>>::set_data(Simplex<1> *abs, Simplex<2> *triang){
	computed = false;
	abscissa = abs;
	triangle = triang;
	clear_all();
};

bool ComputeIntersection<Simplex<1>, Simplex<2>>::is_computed(){
	return computed;
};

void ComputeIntersection<Simplex<1>, Simplex<2>>::set_computed(){
	computed = true;
};

bool ComputeIntersection<Simplex<1>, Simplex<2>>::compute(std::vector<IntersectionPoint<1,2>> &IP12s, bool compute_zeros_plucker_products){

	init_plucker_to_compute();
	computed = true;


    //TODO use reference element to get orientation of sides
	if(((*plucker_products[0] > epsilon) && (*plucker_products[1] < -epsilon) && (*plucker_products[2] > epsilon)) ||
	   ((*plucker_products[0] < -epsilon) && (*plucker_products[1] > epsilon) && (*plucker_products[2] < -epsilon))){
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

		/*xprintf(Msg,"Normalni\n");
		theta.print();
		local_triangle.print();
		xprintf(Msg,"Globale:\n");
		global_triangle.print();*/
		IntersectionPoint<1,2> IP(theta,local_triangle,-1,-1,(*plucker_products[0] > 0 ? 1 : 0));
		//IP12s[0] = IP;
		IP12s.push_back(IP);
		return true;
        //TODO try removing the seocond part of the following condition
	}else if(compute_zeros_plucker_products && (fabs(*plucker_products[0]) <= epsilon || fabs(*plucker_products[1]) <= epsilon || fabs(*plucker_products[2]) <= epsilon)){

		DBGMSG("Intersections - Pathologic case.\n");


		// Pokud je 1 nula = jeden patologický průsečík

		// Pokud jsou 2 nuly => průsečík ve vrcholu

		// Pokud jsou 3 nuly - všechny vypočítat

		for(unsigned int i = 0; i < 3;i++){
			//cout << "PP: " << *plucker_products[i] << endl;
			if(fabs(*plucker_products[i]) <= epsilon){
                // starting point of abscissa
				arma::vec3 A = (*abscissa)[0].getPointCoordinates();
				//arma::vec3 B("9 5 0");
                // direction vector of abscissa
				arma::vec3 U = plucker_coordinates_abscissa[0]->get_u_vector();
				arma::vec3 C = (*triangle)[i][i%2].getPointCoordinates();
				//arma::vec3 C = (*triangle)[i][0].getPointCoordinates();
				//arma::vec3 D("4 4 0");
                // direction vector of triangle side
				arma::vec3 V = plucker_coordinates_triangle[i]->get_u_vector();
				arma::vec3 K = C - A;
                // normal vector to common plane of U and V
				arma::vec3 Det = -arma::cross(U,V);
                
                // we solve following equation for parameters s,t:
                /* A + sU = C + tV
                 * sU - tV = C - A
                 */
                
                //TODO armadillo function for max ??
				//Det.print();
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
				//abscissa is parallel to triangle side
				//TODO compare with epsilon (~ rounding error)
				if(maximum == 0){
					//continue;
					return false;
				}

				double DetX = K[(max_index+2)%3]*V[(max_index+1)%3]
							-K[(max_index+1)%3]*V[(max_index+2)%3];

			    double DetY = K[(max_index+2)%3]*U[(max_index+1)%3]
							-K[(max_index+1)%3]*U[(max_index+2)%3];

				double s;   //parameter on abscissa
				double t;   //parameter on triangle side

				s = DetX/Det[max_index];
				t = DetY/Det[max_index];

                //TODO get from reference element
				if(i == 1){
					t = -t;
				}

				//cout << "s,t: " << s << "," << t << endl;
				//xprintf(Msg, "s,t : %f,%f \n",s,t);

				// IP is outside of triangle side
				if(t > 1+epsilon || t < -epsilon){// || s > 1+epsilon || s < -epsilon){
					//xprintf(Msg,"hoohoo\n");
					//return false;
				}else{

					s = (fabs(s) < epsilon ? 0 : (fabs(1-s) < epsilon ? 1 : s));
					t = (fabs(t) < epsilon ? 0 : (fabs(1-t) < epsilon ? 1 : t));

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

					/*xprintf(Msg,"patologicky\n");
					local_abscissa.print();
					l_triangle.print();*/

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

void ComputeIntersection<Simplex<1>, Simplex<2>>::to_string_plucker_coordinates(){
	cout << "\tPluckerCoordinates Abscissa[0]";
		if(plucker_coordinates_abscissa[0] == NULL){
			cout << "NULL" << endl;
		}else{
			plucker_coordinates_abscissa[0]->toString();
		}
	for(unsigned int i = 0; i < 3;i++){
		cout << "\tPluckerCoordinates Triangle[" << i << "]";
		if(plucker_coordinates_triangle[i] == NULL){
			cout << "NULL" << endl;
		}else{
			plucker_coordinates_triangle[i]->toString();
		}
	}
};

void ComputeIntersection<Simplex<1>, Simplex<2>>::set_plucker_product(double *number, unsigned int i){
	plucker_products[i] = number;
};

double* ComputeIntersection<Simplex<1>, Simplex<2>>::get_plucker_product(unsigned int i){
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

	plucker_coordinates_abscissa.reserve(1);
	plucker_coordinates_tetrahedron.reserve(6);

	plucker_coordinates_abscissa[0] = new Plucker();

	for(unsigned int i = 0; i < 6;i++){
		plucker_coordinates_tetrahedron[i] = new Plucker();
	}

	//plucker_coordinates_tetrahedron.assign(6, new Plucker());
	//plucker_coordinates_abscissa.assign(1, new Plucker());

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

	for(unsigned int j = 0; j < 4;j++){ // for each side of tetrahedron
		for(unsigned int i = 0; i < 3;i++){ // for each side of triangle
			CI12[j].set_pc_triangle(plucker_coordinates_tetrahedron[RefElement<3>::side_lines[j][i]], i);
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
	double epsilon = 64*numeric_limits<double>::epsilon();

    // loop over sides of tetrahedron 
	for(unsigned int side = 0;side < RefElement<3>::n_sides && pocet_pruniku < 2;side++){

        // update plucker products of the side
        //TODO depends on reference element: loop over edges of the side; 
        // possibly can be removed after passing plucker products
		for(unsigned int j = 0; j < side;j++){
			CI12[side].set_plucker_product(CI12[j].get_plucker_product(side-1),j);
		}
		if(!CI12[side].is_computed() // if not computed yet
            && CI12[side].compute(IP12s, true)){    // compute; if intersection exists then continue
			//xprintf(Msg,"Prunik13\n");

			if(IP12s.back().is_patological()){   // resolve pathologic cases
				// Nastavování stěn, které se už nemusí počítat

                //TODO depends on reference element
				if(side == 0){
					if(IP12s.back().get_local_coords2()[0] == 1){
						CI12[1].set_computed();
						CI12[2].set_computed();
					}else if(IP12s.back().get_local_coords2()[1] == 1){
						CI12[1].set_computed();
						CI12[3].set_computed();
					}else if(IP12s.back().get_local_coords2()[2] == 1){
						CI12[2].set_computed();
						CI12[3].set_computed();
					}else{
						CI12[IP12s.back().get_side2() + 1].set_computed();
					}
				}else if(side == 1 && IP12s.back().get_local_coords2()[2] == 1){
					CI12[2].set_computed();
					CI12[3].set_computed();
				}else{
					CI12[IP12s.back().get_side2() + 1].set_computed();
				}
				IP12s.back().set_side2(side);
			}else{
				IP12s.back().set_side2(side);
			}
			pocet_pruniku++;

			//if((IP.get_local_coords1())[0] <= 1 && (IP.get_local_coords1())[0] >= 0){
// 				IP12s.back().print();
				IntersectionPoint<1,3> IP13(IP12s.back());
// 				IP13.print();
				IP13s.push_back(IP13);

			//}
		}
	}

	// Kontrola vytvořených průniků => zda-li není potřeba interpolovat + zda-li se o prunik vubec nejedna:
	if(pocet_pruniku == 1){
		double f_theta = IP13s[IP13s.size()-1].get_local_coords1()[1];
		if(f_theta > 1 || f_theta < 0){
			pocet_pruniku = 0;
			IP13s.pop_back();
		}

	}else if(pocet_pruniku > 1){
        
        /* simplitfy like in NGH
         * intersection.cpp:796
         */
        
		double first_theta = IP13s[IP13s.size()-2].get_local_coords1()[1];
		double second_theta = IP13s[IP13s.size()-1].get_local_coords1()[1];

        //TODO translate comments
        // - compare theta with 1,0 everywhere (without epsilon)
		  // Nejedná se o průnik - celá usečka leží mimo čtyřstěn
		if(((first_theta > 1) && (second_theta > 1)) ||
		   ((first_theta < 0) && (second_theta < 0))){

			pocet_pruniku = 0;
			IP13s.pop_back();
			IP13s.pop_back();
		}else{
            DBGMSG("Intersection interpolation.\n");
			// Jedná se o průnik
			// První souřadnice leží uvnitř čtyřstěnu
			if(first_theta > 1+epsilon || first_theta < -epsilon){
				double theta = first_theta > 1 ? 1 : 0;
				arma::vec::fixed<4> interpolovane = RefElement<3>::line_barycentric_interpolation(IP13s[IP13s.size()-2].get_local_coords2(), IP13s[IP13s.size()-1].get_local_coords2(), first_theta, second_theta,theta);
				arma::vec::fixed<2> inter; inter[0] = 1 - theta; inter[1] = theta;
				IntersectionPoint<1,3> IP13(inter, interpolovane,-1,IP13s[IP13s.size()-2].get_side2(),IP13s[IP13s.size()-2].get_orientation(),true, IP13s[IP13s.size()-2].is_patological());
				IP13s[IP13s.size()-2] = IP13;

				first_theta = theta;
			}else if(fabs(1-first_theta) < epsilon || fabs(first_theta) < epsilon){
				// Hraniční body
				IP13s[IP13s.size()-2].set_is_vertex(true);
				IP13s[IP13s.size()-2].set_is_patological(true);
			}
			// Druhá souřadnice leží uvnitř čtyřstěnu
			if(second_theta > 1+epsilon || second_theta < -epsilon){
				double theta2 = second_theta > 1 ? 1 : 0;
				arma::vec::fixed<2> inter2; inter2[0] = 1 - theta2; inter2[1] = theta2;
				arma::vec::fixed<4> interpolovane2 = RefElement<3>::line_barycentric_interpolation(IP13s[IP13s.size()-2].get_local_coords2(), IP13s[IP13s.size()-1].get_local_coords2(), first_theta, second_theta,theta2);
				IntersectionPoint<1,3> IP13(inter2, interpolovane2,-1,IP13s[IP13s.size()-1].get_side2(),IP13s[IP13s.size()-1].get_orientation(),true, IP13s[IP13s.size()-1].is_patological());
				IP13s[IP13s.size()-1] = IP13;
			}else if(fabs(1-second_theta) < epsilon || fabs(second_theta) < epsilon){
				// hraniční body
				IP13s[IP13s.size()-1].set_is_vertex(true);
				IP13s[IP13s.size()-1].set_is_patological(true);
			}
		}
	}
	return pocet_pruniku;

};

void ComputeIntersection<Simplex<1>, Simplex<3>>::to_string_plucker_coordinates(){
		cout << "\tPluckerCoordinates Abscissa[0]";
		if(plucker_coordinates_abscissa[0] == NULL){
			cout << "NULL" << endl;
		}else{
			plucker_coordinates_abscissa[0]->toString();
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

void ComputeIntersection<Simplex<1>, Simplex<3>>::to_string_plucker_coordinates_tree(){
	cout << "ComputeIntersection<Simplex<1>, <Simplex<3>> Plucker Coordinates Tree:" << endl;
		to_string_plucker_coordinates();
		for(unsigned int i = 0; i < 4;i++){
			cout << "ComputeIntersection<Simplex<1>, Simplex<2>>["<< i <<"] Plucker Coordinates:" << endl;
			CI12[i].to_string_plucker_coordinates();
		}
};

void ComputeIntersection<Simplex<1>, Simplex<3>>::set_plucker_product(double* number, unsigned int index_CI, unsigned index_edge){
	CI12[index_CI].set_plucker_product(number, index_edge);
};

double* ComputeIntersection<Simplex<1>, Simplex<3>>::get_plucker_product(unsigned int index_CI, unsigned index_edge){
	return CI12[index_CI].get_plucker_product(index_edge);
};

/****************************************************************
 * METODY PRO SIMPLEX 2 A SIMPLEX 3
 ****************************************************************/
ComputeIntersection<Simplex<2>, Simplex<3>>::ComputeIntersection(){
	triange = NULL;
	tetrahedron = NULL;

	plucker_coordinates_triangle.reserve(3);
	plucker_coordinates_tetrahedron.reserve(6);
};


ComputeIntersection<Simplex<2>, Simplex<3>>::ComputeIntersection(Simplex<2> &triangle, Simplex<3> &tetr){
	triange = &triangle;
	tetrahedron = &tetr;

	plucker_coordinates_triangle.reserve(3);
	plucker_coordinates_tetrahedron.reserve(6);

    // set CI object for 1D-2D intersection 'tetrahedron edge - triangle'
	for(unsigned int i = 0; i < 6;i++){
		plucker_coordinates_tetrahedron[i] = new Plucker();
		CI12[i].set_data(&tetrahedron->getAbscissa(i), triange);
	}
	// set CI object for 1D-3D intersection 'triangle side - tetrahedron'
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

    // set pointers to Plucker coordinates for 1D-2D
	for(unsigned int i = 0; i < 6;i++){
		CI12[i].set_pc_abscissa(plucker_coordinates_tetrahedron[i]);
		for(unsigned int j = 0; j < 3;j++){
			CI12[i].set_pc_triangle(plucker_coordinates_triangle[j],j);
		}
	}
	
	// set pointers to Plucker coordinates for 1D-3D
	for(unsigned int i = 0; i < 3;i++){
		CI13[i].set_pc_abscissa(plucker_coordinates_triangle[i]);
		for(unsigned int j = 0; j < 6;j++){
			CI13[i].set_pc_tetrahedron(plucker_coordinates_tetrahedron[j],j);
		}
		CI13[i].init();
	}
};

void ComputeIntersection<Simplex<2>, Simplex<3>>::compute(IntersectionPolygon &lokalni_mnohouhelnik){

	std::vector<IntersectionPoint<1,2>> IP12s;
	std::vector<IntersectionPoint<1,3>> IP13s;
	int pocet_pruniku = 0;
	int pocet_13_pruniku;

	//cout << "ComputeIntersection<Simplex<2>, Simplex<3>>::compute - edges triangle vs tetrahedron" << endl;
		for(unsigned int i = 0; i < 3;i++){
			pocet_13_pruniku = CI13[(3-i)%3].compute(IP13s);
//             DBGMSG("CI23: number of 1-3 intersections = %d\n",pocet_13_pruniku);
			//(triange->getAbscissa(i)).toString();
			// Vždy by měl být počet průniku 2 nebo 0
			if(pocet_13_pruniku == 1){
				IP13s[IP13s.size() - 1].set_side1((3-i)%3);
				IntersectionPoint<3,1> IP31(IP13s[IP13s.size() - 1]);
				IntersectionPoint<3,2> IP32(IP31);
				IntersectionPoint<2,3> IP23(IP32);
				lokalni_mnohouhelnik.add_ipoint(IP23);

			}else if(pocet_13_pruniku == 2){
//                 DBGMSG("i=%d,  set side1 %d, side2 %d\n",i, (3-i)%3, IP13s[IP13s.size() - 2].get_side2());
				IP13s[IP13s.size() - 2].set_side1((3-i)%3);
				//IntersectionPoint<3,1> IP31 = IntersectionLocal::flipDimension<3,1>(IP13s[IP13s.size() - 2]);
				IntersectionPoint<3,1> IP31(IP13s[IP13s.size() - 2]);
				//IntersectionPoint<3,2> IP32 = IntersectionLocal::interpolateDimension<3,2>(IP31);
				IntersectionPoint<3,2> IP32(IP31);
				IntersectionPoint<2,3> IP23(IP32);
				lokalni_mnohouhelnik.add_ipoint(IP23);

//                 DBGMSG("i=%d,  set side1 %d, side2 %d\n",i, (3-i)%3, IP13s[IP13s.size() - 1].get_side2());
				IP13s[IP13s.size() - 1].set_side1((3-i)%3);
				IntersectionPoint<3,1> IP31_2(IP13s[IP13s.size() - 1]);
				IntersectionPoint<3,2> IP32_2(IP31_2);
				IntersectionPoint<2,3> IP23_2(IP32_2);
				lokalni_mnohouhelnik.add_ipoint(IP23_2);

			}else if(pocet_13_pruniku > 2){
				// TODO: nahradit Assertem
				xprintf(Msg, "JINY POCET PRUNIKU\n");
			}


		}
	// Optimalizace: znovu použití již vypočítaných součinů


	for(unsigned int i = 0; i < 3;i++){
		CI12[0].set_plucker_product(CI13[i].get_plucker_product(0,0),i);
		CI12[1].set_plucker_product(CI13[i].get_plucker_product(0,1),i);
		CI12[2].set_plucker_product(CI13[i].get_plucker_product(0,2),i);
		CI12[3].set_plucker_product(CI13[i].get_plucker_product(1,1),i);
		CI12[4].set_plucker_product(CI13[i].get_plucker_product(1,2),i);
		CI12[5].set_plucker_product(CI13[i].get_plucker_product(2,2),i);
	}

	//cout << "ComputeIntersection<Simplex<2>, Simplex<3>>::compute - edges tetrahedron vs triangle" << endl;
	//double epsilon = 64*numeric_limits<double>::epsilon();
	for(unsigned int i = 0; i < 6;i++){
			if(CI12[i].compute(IP12s, false)){
				pocet_pruniku++;
				IP12s.back().set_side1(i);
				if((IP12s.back().get_local_coords1())[0] <= 1 && (IP12s.back().get_local_coords1())[0] >= 0){
					/*if(fabs(1-IP12s.back().get_local_coords1()[0]) < epsilon || fabs(IP12s.back().get_local_coords1()[0]) < epsilon){
						IP12s.back().setIsPatological(true);
					}*/
					//IP.print();
// 					DBGMSG("CI23: number of 1-2 intersections = %d\n",pocet_pruniku);
//                     DBGMSG("i=%d,  set side1 %d, side2 %d ori %d\n",i, IP12s.back().get_side1(), IP12s.back().get_side2(), IP12s.back().get_orientation());
									IntersectionPoint<2,1> IP21(IP12s.back());
									IntersectionPoint<2,3> IP23(IP21);
									//IntersectionPoint<2,3> IP23 = IntersectionLocal::interpolateDimension<2,2>(IP22);
									//IP23.print();
									lokalni_mnohouhelnik.add_ipoint(IP23);


				}

				//IP12s.pop_back();
			}
		}


};

void ComputeIntersection<Simplex<2>, Simplex<3>>::to_string_plucker_coordinates(){
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

void ComputeIntersection<Simplex<2>, Simplex<3>>::to_string_plucker_coordinates_tree(){
	cout << "ComputeIntersection<Simplex<2>, <Simplex<3>> Plucker Coordinates Tree:" << endl;
	to_string_plucker_coordinates();
	for(unsigned int i = 0; i < 6;i++){
		cout << "ComputeIntersection<Simplex<1>, Simplex<2>>["<< i <<"] Plucker Coordinates:" << endl;
		CI12[i].to_string_plucker_coordinates();
	}
	for(unsigned int i = 0; i < 3;i++){
		CI13[i].to_string_plucker_coordinates_tree();
	}
};

} // END namespace_close
