/*
 *      Author: viktor
 */

#include "computeintersection.h"
#include "mesh/ref_element.hh"
#include "system/system.hh"

#include "plucker.h"
#include "intersectionpoint.h"
#include "intersectionpolygon.h"

using namespace std;
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
		plucker_coordinates_abscissa[0]->compute((*abscissa)[0].point_coordinates(),
												 (*abscissa)[1].point_coordinates());
	}
	// if not already computed, compute plucker coordinates of triangle sides
	for(unsigned int side = 0; side < RefElement<2>::n_sides; side++){
		if(!plucker_coordinates_triangle[side]->is_computed()){
			plucker_coordinates_triangle[side]->compute((*triangle)[side][0].point_coordinates(), (*triangle)[side][1].point_coordinates());
		}
	}

// 	DBGMSG("Abscissa:\n");
//     (*abscissa)[0].point_coordinates().print();
//     (*abscissa)[1].point_coordinates().print();
    
	// compute Plucker products (abscissa X triangle side)
	for(unsigned int side = 0; side < RefElement<2>::n_sides; side++){
		if(plucker_products[side] == NULL){
//             (*plucker_coordinates_abscissa[0]).toString();
//             (*plucker_coordinates_triangle[side]).toString();
			plucker_products[side] = new double((*plucker_coordinates_abscissa[0])*(*plucker_coordinates_triangle[side]));
            
//             DBGMSG("triangle side:\n");
//             (*triangle)[side][0].point_coordinates().print();
//             (*triangle)[side][1].point_coordinates().print();
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
		arma::vec3 vec = (*abscissa)[1].point_coordinates() - (*abscissa)[0].point_coordinates();
		unsigned int i = 0;
		double max = vec[0];

		if(fabs((double)vec[1]) > fabs(max)){ max = vec[1]; i = 1;}
		if(fabs((double)vec[2]) > fabs(max)){ max = vec[2]; i = 2;}

		arma::vec3 global_triangle =
		local_triangle[0]*(*triangle)[0][0].point_coordinates() +
		local_triangle[1]*(*triangle)[0][1].point_coordinates() +
		local_triangle[2]*(*triangle)[1][1].point_coordinates();
		theta[1] = (-(*abscissa)[0].point_coordinates()[i] + global_triangle[i])/max;
		theta[0] = 1 - theta[1];

		/*xprintf(Msg,"Normalni\n");
		theta.print();
		local_triangle.print();
		xprintf(Msg,"Globale:\n");
		global_triangle.print();*/
//         IntersectionPoint<1,2> IP(theta,local_triangle,-1,-1,(*plucker_products[0] > 0 ? 1 : 0));
        /*TODO: we do not test whether the IP is not at the vertex of abscissa:
         * we solve intersection of triangle and line (not abscissa)
         * this should be done (after this function) when solving 1d-2d
         */
		IntersectionPoint<1,2> IP(theta,local_triangle);
        IP.set_topology(0, 1, 0, 2);
        IP.set_plucker_flags((*plucker_products[0] > 0 ? 1 : 0), false);
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
				arma::vec3 A = (*abscissa)[0].point_coordinates();
				//arma::vec3 B("9 5 0");
                // direction vector of abscissa
				arma::vec3 U = plucker_coordinates_abscissa[0]->get_u_vector();
				arma::vec3 C = (*triangle)[i][i%2].point_coordinates();
				//arma::vec3 C = (*triangle)[i][0].point_coordinates();
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

				// IP is outside of triangle side
				if(t > 1+epsilon || t < -epsilon){// || s > 1+epsilon || s < -epsilon){
					//xprintf(Msg,"hoohoo\n");
					//return false;
				}else{
                    IntersectionPoint<1,2> IP;
                    IP.set_plucker_flags(1, true);  // orientation and pathologic flag
                    
                    // possibly set abscissa vertex {0,1}
                    if( fabs(s) < epsilon)       { s = 0; IP.set_topology_A(0,0);}
                    else if(fabs(1-s) < epsilon) { s = 1; IP.set_topology_A(1,0);}
                    else                         {        IP.set_topology_A(0,1);}   // no vertex, line 0, dim = 1
                    
                    // possibly set triangle vertex {0,1,2}
                    if( fabs(t) < epsilon)       { t = 0; IP.set_topology_B(RefElement<2>::interact<0,1>(i)[0],0);}
                    else if(fabs(1-t) < epsilon) { t = 1; IP.set_topology_B(RefElement<2>::interact<0,1>(i)[1],0);}
                    else                         {        IP.set_topology_B(i,1);}   // no vertex, side i, dim = 1
                    
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
					 * TODO: according to ref element
					 * */
					l_triangle[(3-i)%3] = 1 - t;
					l_triangle[(4-i)%3] = t;
					l_triangle[2-i] = 0;

// 					local_abscissa.print();
// 					l_triangle.print();
                    
					IP.set_coordinates(local_abscissa,l_triangle);
                    
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

void ComputeIntersection<Simplex<1>, Simplex<2>>::print_plucker_coordinates(std::ostream &os){
	os << "\tPluckerCoordinates Abscissa[0]";
		if(plucker_coordinates_abscissa[0] == NULL){
			os << "NULL" << endl;
		}else{
			os << plucker_coordinates_abscissa[0];
		}
	for(unsigned int i = 0; i < 3;i++){
		os << "\tPluckerCoordinates Triangle[" << i << "]";
		if(plucker_coordinates_triangle[i] == NULL){
			os << "NULL" << endl;
		}else{
			os << plucker_coordinates_triangle[i];
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
			CI12[j].set_pc_triangle(plucker_coordinates_tetrahedron[RefElement<3>::interact<1,2>(j)[i]], i); //[RefElement<3>::side_lines[j][i]], i);
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

unsigned int ComputeIntersection<Simplex<1>, Simplex<3>>::compute(std::vector<IntersectionPoint<1,3>> &IP13s){

	std::vector<IntersectionPoint<1,2>> IP12s;
	unsigned int pocet_pruniku = 0;

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

            IntersectionPoint<1,2> IP = IP12s.back();   // shortcut
            IntersectionPoint<1,3> IP13(IP, side);
        
			if(IP.is_pathologic()){   // resolve pathologic cases
                
                // set the 'computed' flag on the connected sides by IP
                if(IP.dim_B() == 0) // IP is vertex of triangle
                {
                    for(unsigned int s=0; s < 3; s++)   //sides per node
                        CI12[RefElement<3>::interact<2,0>(IP.idx_B())[s]].set_computed();
                    IP13.set_topology_B(RefElement<3>::interact<0,2>(side)[IP.idx_B()], IP.dim_B());
                }
                if(IP.dim_B() == 1) // IP is on edge of triangle
                {
                    for(unsigned int s=0; s < RefElement<3>::n_sides_per_line; s++)
                        CI12[RefElement<3>::interact<2,1>(IP.idx_B())[s]].set_computed();
                    IP13.set_topology_B(RefElement<3>::interact<1,2>(side)[IP.idx_B()], IP.dim_B());
                }
			}
			
			pocet_pruniku++;
            IP13s.push_back(IP13);
		}
	}
    
    ASSERT_LESS(pocet_pruniku,3);
    
	// in the case, that line goes through vertex, but outside tetrahedron (touching vertex)
	if(pocet_pruniku == 1){
		double f_theta = IP13s[IP13s.size()-1].local_bcoords_A()[1];
		if(f_theta > 1 || f_theta < 0){
			pocet_pruniku = 0;
			IP13s.pop_back();
		}

	}else if(pocet_pruniku == 2){
        
        // swap the ips in the line coordinate direction (0->1 : a1->a2)
        IntersectionPoint<1,3> *a1, *a2;
        double t1, t2;
        if(IP13s[IP13s.size()-2].local_bcoords_A()[1] > IP13s[IP13s.size()-1].local_bcoords_A()[1])
        {
            a1 = &(IP13s[IP13s.size()-1]);
            a2 = &(IP13s[IP13s.size()-2]);
        }
        else
        {
            a1 = &(IP13s[IP13s.size()-2]);
            a2 = &(IP13s[IP13s.size()-1]);
        }
        // get first and second theta (coordinate of ips on the line)
        t1 = a1->local_bcoords_A()[1];
        t2 = a2->local_bcoords_A()[1];
        
        // cut off the line by the abscissa points
        if(t1 < 0) t1 = 0;
        if(t2 > 1) t2 = 1;
        
        if(t2 < t1) { // then the intersection is outside the abscissa => NO intersection
            pocet_pruniku = 0;
            IP13s.pop_back();
            IP13s.pop_back(); 
            return pocet_pruniku;
        }
        
        if(t1 == 0) // interpolate IP a1
        {
            arma::vec::fixed<4> interpolovane = RefElement<3>::line_barycentric_interpolation(a1->local_bcoords_B(), 
                                                                                              a2->local_bcoords_B(), 
                                                                                              a1->local_bcoords_A()[1],
                                                                                              a2->local_bcoords_A()[1], 
                                                                                              t1);
            arma::vec::fixed<2> inter({1 - t1, t1});    // barycentric coords
            a1->set_coordinates(inter,interpolovane);
//             a1->set_topology_EE(unset_loc_idx, a1->is_pathologic());  // edge index is set later
            // here we can set only local index of the vertex on the line
            DBGMSG("E-E 0\n");
            a1->set_topology(0, 0, 0,3);
        }
        if(t2 == 1) // interpolate IP a2
        {
            arma::vec::fixed<4> interpolovane = RefElement<3>::line_barycentric_interpolation(a1->local_bcoords_B(), 
                                                                                              a2->local_bcoords_B(), 
                                                                                              a1->local_bcoords_A()[1],
                                                                                              a2->local_bcoords_A()[1], 
                                                                                              t2);
            arma::vec::fixed<2> inter({1 - t2, t2});      // barycentric coords
            a2->set_coordinates(inter,interpolovane);
//             a2->set_topology_EE(unset_loc_idx, a2->is_pathologic());  // edge index is set later
            // here we can set only local index of the vertex on the line
            DBGMSG("E-E 1\n");
            a2->set_topology(1, 0, 0,3);
        }
    }
    return pocet_pruniku;
};

void ComputeIntersection<Simplex<1>, Simplex<3>>::print_plucker_coordinates(std::ostream &os){
		os << "\tPluckerCoordinates Abscissa[0]";
		if(plucker_coordinates_abscissa[0] == NULL){
			os << "NULL" << endl;
		}else{
			os << plucker_coordinates_abscissa[0];
		}

	for(unsigned int i = 0; i < 6;i++){
		os << "\tPluckerCoordinates Tetrahedron[" << i << "]";
		if(plucker_coordinates_tetrahedron[i] == NULL){
			os << "NULL" << endl;
		}else{
			os << plucker_coordinates_tetrahedron[i];
		}
	}
};

void ComputeIntersection<Simplex<1>, Simplex<3>>::print_plucker_coordinates_tree(std::ostream &os){
	os << "ComputeIntersection<Simplex<1>, <Simplex<3>> Plucker Coordinates Tree:" << endl;
		print_plucker_coordinates(os);
		for(unsigned int i = 0; i < 4;i++){
			os << "ComputeIntersection<Simplex<1>, Simplex<2>>["<< i <<"] Plucker Coordinates:" << endl;
			CI12[i].print_plucker_coordinates(os);
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
		CI12[i].set_data(&tetrahedron->abscissa(i), triange);
	}
	// set CI object for 1D-3D intersection 'triangle side - tetrahedron'
	for(unsigned int i = 0; i < 3;i++){
		plucker_coordinates_triangle[i] = new Plucker();
		CI13[i].set_data(&triange->abscissa(i) , tetrahedron);
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
	unsigned int pocet_13_pruniku;

	for(unsigned int i = 0; i < RefElement<2>::n_lines;i++){    // go through triangle lines
		pocet_13_pruniku = CI13[(3-i)%3].compute(IP13s);
        ASSERT(pocet_13_pruniku < 3, "Impossible number of intersection.");
//         DBGMSG("CI23: number of 1-3 intersections = %d\n",pocet_13_pruniku);
        
        for(unsigned int n=1; n <= pocet_13_pruniku; n++){
            IntersectionPoint<1,3> IP (IP13s[IP13s.size()-n]);
            
            IntersectionPoint<3,1> IP31(IP);             // switch idx_A and idx_B and coords
            IntersectionPoint<3,2> IP32(IP31, (3-i)%3);  // interpolation uses local_bcoords_B and given idx_B
            IntersectionPoint<2,3> IP23(IP32);           // switch idx_A and idx_B and coords back
            
            if( IP.dim_A() == 0 ) // if IP is vertex of triangle
            {
                // we are on line (3-i)%3 of the triangle, and IP.idx_A contains local node of the line
                // E-E, we know vertex index
                IP23.set_topology_A(RefElement<2>::interact<0,1>((3-i)%3)[IP.idx_A()], 0);
            }
            
            lokalni_mnohouhelnik.add_ipoint(IP23);
        }
    }

    // Optimalization: reusage of the Plucker coordinates products already computed
	for(unsigned int i = 0; i < 3;i++){
		CI12[0].set_plucker_product(CI13[i].get_plucker_product(0,0),i);
		CI12[1].set_plucker_product(CI13[i].get_plucker_product(0,1),i);
		CI12[2].set_plucker_product(CI13[i].get_plucker_product(0,2),i);
		CI12[3].set_plucker_product(CI13[i].get_plucker_product(1,1),i);
		CI12[4].set_plucker_product(CI13[i].get_plucker_product(1,2),i);
		CI12[5].set_plucker_product(CI13[i].get_plucker_product(2,2),i);
	}

	for(unsigned int i = 0; i < 6;i++){
		if(CI12[i].compute(IP12s, false)){
            IntersectionPoint<1,2> IP = IP12s.back();
            // Check whether the IP is on the abscissa (side of triangle)
			if((IP.local_bcoords_A())[0] <= 1 && (IP.local_bcoords_A())[0] >= 0){
// 				DBGMSG("CI23: number of 1-2 intersections = %d\n",pocet_pruniku);
				IntersectionPoint<2,1> IP21(IP);
				IntersectionPoint<2,3> IP23(IP21,i);
                
                if(IP.dim_A() == 0) // IP is vertex of line (i.e. line of tetrahedron)
                    IP23.set_topology_B(RefElement<3>::interact<0,1>(i)[IP.idx_A()], 0);
                
				//IP23.print();
				lokalni_mnohouhelnik.add_ipoint(IP23);
			}
		}
	}
};

void ComputeIntersection<Simplex<2>, Simplex<3>>::print_plucker_coordinates(std::ostream &os){
	for(unsigned int i = 0; i < 3;i++){
		os << "\tPluckerCoordinates Triangle[" << i << "]";
		if(plucker_coordinates_triangle[i] == NULL){
			os << "NULL" << endl;
		}else{
			os << plucker_coordinates_triangle[i];
		}
	}
	for(unsigned int i = 0; i < 6;i++){
		os << "\tPluckerCoordinates Tetrahedron[" << i << "]";
		if(plucker_coordinates_tetrahedron[i] == NULL){
			os << "NULL" << endl;
		}else{
			os << plucker_coordinates_tetrahedron[i];
		}
	}
};

void ComputeIntersection<Simplex<2>, Simplex<3>>::print_plucker_coordinates_tree(std::ostream &os){
	os << "ComputeIntersection<Simplex<2>, <Simplex<3>> Plucker Coordinates Tree:" << endl;
	print_plucker_coordinates(os);
	for(unsigned int i = 0; i < 6;i++){
		os << "ComputeIntersection<Simplex<1>, Simplex<2>>["<< i <<"] Plucker Coordinates:" << endl;
		CI12[i].print_plucker_coordinates(os);
	}
	for(unsigned int i = 0; i < 3;i++){
		CI13[i].print_plucker_coordinates_tree(os);
	}
};

} // END namespace_close
