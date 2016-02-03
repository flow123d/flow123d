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

/*************************************************************************************************************
 *                                  COMPUTE INTERSECTION FOR:             1D AND 2D
 ************************************************************************************************************/
ComputeIntersection<Simplex<1>, Simplex<2>>::ComputeIntersection()
: computed_(false), abscissa_(nullptr), triangle_(nullptr)
{
    plucker_coordinates_abscissa_ = nullptr;
	plucker_coordinates_triangle_.resize(3);
};

ComputeIntersection<Simplex<1>, Simplex<2>>::ComputeIntersection(computeintersection::Simplex< 1 >& abscissa, 
                                                                 computeintersection::Simplex< 2 >& triangle)
: computed_(false), abscissa_(&abscissa), triangle_(&triangle)
{
    plucker_coordinates_abscissa_ = new Plucker();
    plucker_coordinates_triangle_.resize(3);

	for(unsigned int i = 0; i < 3;i++){
		plucker_coordinates_triangle_[i] = new Plucker();
	}

	clear_all();
};

void ComputeIntersection<Simplex<1>, Simplex<2>>::clear_all(){
	for(unsigned int i = 0; i < 3;i++)
		plucker_products_[i] = nullptr;
};

void ComputeIntersection<Simplex<1>, Simplex<2>>::compute_plucker_products(){

    // if not already computed, compute plucker coordinates of abscissa
	if(!plucker_coordinates_abscissa_->is_computed()){
		plucker_coordinates_abscissa_->compute((*abscissa_)[0].point_coordinates(),
											  (*abscissa_)[1].point_coordinates());
	}
	// if not already computed, compute plucker coordinates of triangle sides
	for(unsigned int side = 0; side < RefElement<2>::n_sides; side++){
		if(!plucker_coordinates_triangle_[side]->is_computed()){
			plucker_coordinates_triangle_[side]->compute((*triangle_)[side][0].point_coordinates(), 
                                                        (*triangle_)[side][1].point_coordinates());
		}
	}

// 	DBGMSG("Abscissa:\n");
//     (*abscissa_)[0].point_coordinates().print();
//     (*abscissa_)[1].point_coordinates().print();
    
	// compute Plucker products (abscissa X triangle side)
	for(unsigned int side = 0; side < RefElement<2>::n_sides; side++){
		if(plucker_products_[side] == nullptr){
//             (*plucker_coordinates_abscissa_).toString();
//             (*plucker_coordinates_triangle_[side]).toString();
			plucker_products_[side] = new double((*plucker_coordinates_abscissa_)*(*plucker_coordinates_triangle_[side]));
            
//             DBGMSG("triangle side:\n");
//             (*triangle_)[side][0].point_coordinates().print();
//             (*triangle_)[side][1].point_coordinates().print();
		}
// 		DBGMSG("Plucker product = %f\n", *(plucker_products_[side]));
	}

};

void ComputeIntersection<Simplex<1>, Simplex<2>>::set_data(computeintersection::Simplex< 1 >& abscissa, 
                                                           computeintersection::Simplex< 2 >& triangle){
	computed_ = false;
	abscissa_ = &abscissa;
	triangle_ = &triangle;
	clear_all();
};

bool ComputeIntersection< Simplex< 1  >, Simplex< 2  > >::compute_plucker(IntersectionPoint< 1, 2 > &IP)
{
    // compute local barycentric coordinates of IP: see formula (3) on pg. 12 in BP VF
    // local alfa = w2/sum; local beta = w1/sum; => local barycentric coordinates in the triangle
    // where sum = w0+w1+w2
    
    double w_sum = 0;
    // plucker products with correct signs according to ref_element, ordered by sides
    double w[3] = {signed_plucker_product(0),
                   signed_plucker_product(1),
                   signed_plucker_product(2)};
        
    // local barycentric coordinates of IP, depends on barycentric coordinates order !!!
    arma::vec3 local_triangle({w[2],w[1],w[0]});
    local_triangle = local_triangle / (w[0]+w[1]+w[2]);
    
    // local coordinate T on the line
    // for i-th coordinate it holds: (from formula (4) on pg. 12 in BP VF)
    // T = localAbscissa= (- A(i) + ( 1 - alfa - beta ) * V0(i) + alfa * V1(i) + beta * V2 (i)) / U(i)
    // let's choose [max,i] = max {U(i)}
    arma::vec3 u = (*abscissa_)[1].point_coordinates() - (*abscissa_)[0].point_coordinates();
    
    //find max in u in abs value:
    unsigned int i = 0; // index of maximum in u
    double max = u[0];  // maximum in u
    if(fabs((double)u[1]) > fabs(max)){ max = u[1]; i = 1;}
    if(fabs((double)u[2]) > fabs(max)){ max = u[2]; i = 2;}

    // global coordinates in triangle
    arma::vec3 global_triangle =
    local_triangle[0]*(*triangle_)[0][0].point_coordinates() +
    local_triangle[1]*(*triangle_)[0][1].point_coordinates() +
    local_triangle[2]*(*triangle_)[1][1].point_coordinates();

    //theta on abscissa
    double t =  (-(*abscissa_)[0].point_coordinates()[i] + global_triangle[i])/max;

    /*
    DBGMSG("Coordinates: line; local and global triangle\n");
    theta.print();
    local_triangle.print();
    global_triangle.print();
    */
    IP.set_topology(0, 1, 0, 2);
    IP.set_orientation(signed_plucker_product(0) > 0 ? 1 : 0);
    arma::vec2 local_abscissa = {1-t, t};

    IP.set_coordinates(local_abscissa,local_triangle);
    return true;
}

bool ComputeIntersection< Simplex< 1  >, Simplex< 2  > >::compute_pathologic(unsigned int side, IntersectionPoint< 1, 2 > &IP)
{
//      DBGMSG("PluckerProduct[%d]: %f\n",side, *plucker_products_[side]);
        if( std::abs(*plucker_products_[side]) <= rounding_epsilon ){// fabs(*plucker_products_[side]) <= rounding_epsilon){
            
            // We solve following equation for parameters s,t:
            /* A + sU = C + tV = intersection point
             * sU - tV = C - A
             * 
             * which is by components:
             * (u1  -v1) (s) = (c1-a1)
             * (u2  -v2) (t) = (c2-a2)
             * (u3  -v3)     = (c3-a3)
             * 
             * these are 3 equations for variables s,t
             * see (4.3) on pg. 19 in DP VF
             * 
             * We will solve this using Crammer's rule for the maximal subdeterminant det_ij of the matrix.
             * s = detX_ij / det_ij
             * t = detY_ij / det_ij
             */
            
            // starting point of abscissa
            arma::vec3 A = (*abscissa_)[0].point_coordinates();
            // direction vector of abscissa
            arma::vec3 U = plucker_coordinates_abscissa_->get_u_vector();
            // vertex of triangle side
            arma::vec3 C = (*triangle_)[side][side%2].point_coordinates();
            // direction vector of triangle side
            arma::vec3 V = plucker_coordinates_triangle_[side]->get_u_vector();
            // right hand side
            arma::vec3 K = C - A;
            // subdeterminants det_ij of the system equal minus normal vector to common plane of U and V
            // det_12 =  (-UxV)[1]
            // det_13 = -(-UxV)[2]
            // det_23 =  (-UxV)[3]
            arma::vec3 Det = -arma::cross(U,V);
            
//             A.print();
//             U.print();
//             C.print();
//             V.print();
//             K.print();
//             Det.print();
//             cout <<endl;

            unsigned int max_index = 0;
            double maximum = fabs(Det[0]);
            if(fabs((double)Det[1]) > maximum){
                maximum = fabs(Det[1]);
                max_index = 1;
            }
            if(fabs((double)Det[2]) > maximum){
                maximum = fabs(Det[2]);
                max_index = 2;
            }
            //abscissa is parallel to triangle side
            if(std::abs(maximum) <= rounding_epsilon) return false;

            // map maximum index in {-UxV} to i,j of subdeterminants
            //              i j
            // max_index 0: 1 2
            //           1: 2 0  (switch  due to sign change)
            //           2: 0 1
            unsigned int i = (max_index+1)%3,
                         j = (max_index+2)%3;
                         
            double DetX = -K[i]*V[j] + K[j]*V[i];
            double DetY = -K[i]*U[j] + K[j]*U[i];

            double s = DetX/Det[max_index];   //parameter on abscissa
            double t = DetY/Det[max_index];   //parameter on triangle side

            // change sign according to side orientation
            if(RefElement<2>::normal_orientation(side)) t=-t;

            //DBGMSG("s = %f; t = %f\n",s,t);

            // if IP is inside of triangle side
            if(t >= -geometry_epsilon && t <= 1+geometry_epsilon){

                IP.set_orientation(2);  // set orientation as a pathologic case ( > 1)

                // possibly set abscissa vertex {0,1}
                if( fabs(s) <= geometry_epsilon)       { s = 0; IP.set_topology_A(0,0);}
                else if(fabs(1-s) <= geometry_epsilon) { s = 1; IP.set_topology_A(1,0);}
                else                         {        IP.set_topology_A(0,1);}   // no vertex, line 0, dim = 1
                
                // possibly set triangle vertex {0,1,2}
                if( fabs(t) <= geometry_epsilon)       { t = 0; IP.set_topology_B(RefElement<2>::interact<0,1>(side)[RefElement<2>::normal_orientation(side)],0);  }
                else if(fabs(1-t) <= geometry_epsilon) { t = 1; IP.set_topology_B(RefElement<2>::interact<0,1>(side)[1-RefElement<2>::normal_orientation(side)],0);}
                else                         {        IP.set_topology_B(side,1);}   // no vertex, side side, dim = 1
                
                arma::vec2 local_abscissa({1-s, s});
                arma::vec3 local_triangle({0,0,0});

                // set local triangle barycentric coordinates according to nodes of the triangle side:
                local_triangle[RefElement<2>::interact<0,1>(side)[RefElement<2>::normal_orientation(side)]] = 1 - t;
                local_triangle[RefElement<2>::interact<0,1>(side)[1-RefElement<2>::normal_orientation(side)]] = t;
//              local_abscissa.print();
//              local_triangle.print();
                
                IP.set_coordinates(local_abscissa,local_triangle);

                return true; // IP found
            }
        }
    return false;   // IP NOT found
}


bool ComputeIntersection<Simplex<1>, Simplex<2>>::compute(std::vector<IntersectionPoint<1,2>> &IP12s, 
                                                          bool compute_zeros_plucker_products){

    compute_plucker_products();
    computed_ = true;
    
    // test whether all plucker products have the same sign
    if(((signed_plucker_product(0) > rounding_epsilon) && (signed_plucker_product(1) > rounding_epsilon) && (signed_plucker_product(2) > rounding_epsilon)) ||
       ((signed_plucker_product(0) < -rounding_epsilon) && (signed_plucker_product(1) < -rounding_epsilon) && (signed_plucker_product(2) < -rounding_epsilon))){
        
        IntersectionPoint<1,2> IP;
        
        if(compute_plucker(IP))
        {
            IP12s.push_back(IP);
            return true;
        }
    }else if(compute_zeros_plucker_products){

        DBGMSG("Intersections - Pathologic case.\n");
        // 1 zero product -> IP is on the triangle side
        // 2 zero products -> IP is at the vertex of triangle (there is no other IP)
        // 3 zero products: 
        //      -> IP is at the vertex of triangle but the line is parallel to opossite triangle side
        //      -> triangle side is part of the line (and otherwise)
        IntersectionPoint<1,2> IP;
        for(unsigned int i = 0; i < 3;i++){
            if(compute_pathologic(i,IP))
            {
                IP12s.push_back(IP);
                return true;
            }
        }
    }
    
    return false;   // if IP not found before
};

bool ComputeIntersection< Simplex< 1  >, Simplex< 2  > >::compute_final(vector< IntersectionPoint< 1, 2 > >& IP12s)
{
    compute_plucker_products();
    computed_ = true;
    
    // test whether all plucker products have the same sign
    if(((signed_plucker_product(0) > rounding_epsilon) && (signed_plucker_product(1) > rounding_epsilon) && (signed_plucker_product(2) > rounding_epsilon)) ||
       ((signed_plucker_product(0) < -rounding_epsilon) && (signed_plucker_product(1) < -rounding_epsilon) && (signed_plucker_product(2) < -rounding_epsilon)))
    {
        IntersectionPoint<1,2> IP;
        compute_plucker(IP);
        
        // if computing 1d-2d intersection as final product, cut the line
        arma::vec::fixed<2> theta;
        double t = IP.local_bcoords_A()[1];
        if(t >= -geometry_epsilon && t <= 1+geometry_epsilon){
                // possibly set abscissa vertex {0,1}
                if( fabs(t) <= geometry_epsilon)       { theta = {1,0}; IP.set_topology_A(0,0);}
                else if(fabs(1-t) <= geometry_epsilon) { theta = {0,1}; IP.set_topology_A(1,0);}
                IP12s.push_back(IP);
                return true;
        }
        else return false;   
    }
    else
    {
        DBGMSG("Intersections - Pathologic case.\n");
        // 1 zero product -> IP is on the triangle side
        // 2 zero products -> IP is at the vertex of triangle (there is no other IP)
        // 3 zero products: 
        //      -> IP is at the vertex of triangle but the line is parallel to opossite triangle side
        //      -> triangle side is part of the line (and otherwise)      
        unsigned int n_found = 0;
        
        for(unsigned int i = 0; i < 3;i++){
            IntersectionPoint<1,2> IP;
            if (compute_pathologic(i,IP))
            {
                double t = IP.local_bcoords_A()[1];
                if(t >= -geometry_epsilon && t <= 1+geometry_epsilon)
                {
                    if(n_found > 0)
                    {
                        // if the IP has been found already
                        if(IP12s.back().local_bcoords_A()[1] == t)
                            continue;
                        
                        // sort the IPs in the direction of the abscissa
                        if(IP12s.back().local_bcoords_A()[1] > t)
                            std::swap(IP12s.back(),IP);
                    }
                    
                    IP12s.push_back(IP);
                    
                    n_found++;
                }
            }
        }
        return n_found > 0;
    }
}


void ComputeIntersection<Simplex<1>, Simplex<2>>::print_plucker_coordinates(std::ostream &os){
	os << "\tPluckerCoordinates Abscissa[0]";
		if(plucker_coordinates_abscissa_ == nullptr){
			os << "NULL" << endl;
		}else{
			os << plucker_coordinates_abscissa_;
		}
	for(unsigned int i = 0; i < 3;i++){
		os << "\tPluckerCoordinates Triangle[" << i << "]";
		if(plucker_coordinates_triangle_[i] == nullptr){
			os << "NULL" << endl;
		}else{
			os << plucker_coordinates_triangle_[i];
		}
	}
};





/*************************************************************************************************************
 *                                  COMPUTE INTERSECTION FOR:             1D AND 3D
 ************************************************************************************************************/
ComputeIntersection<Simplex<1>, Simplex<3>>::ComputeIntersection()
{
    plucker_coordinates_abscissa_ = nullptr;
    plucker_coordinates_tetrahedron.resize(6, nullptr);
};

ComputeIntersection<Simplex<1>, Simplex<3>>::ComputeIntersection(computeintersection::Simplex< 1 >& abscissa, 
                                                                 computeintersection::Simplex< 3 >& tetrahedron)
{
    plucker_coordinates_abscissa_ = new Plucker();
    plucker_coordinates_tetrahedron.resize(6);

	for(unsigned int i = 0; i < 6;i++){
		plucker_coordinates_tetrahedron[i] = new Plucker();
	}

    set_data(abscissa, tetrahedron);
};



void ComputeIntersection<Simplex<1>, Simplex<3>>::init(){
    
	for(unsigned int i = 0; i < 4;i++){
		CI12[i].set_pc_abscissa(plucker_coordinates_abscissa_);
	}

	for(unsigned int j = 0; j < 4;j++){ // for each side of tetrahedron
		for(unsigned int i = 0; i < 3;i++){ // for each side of triangle
			CI12[j].set_pc_triangle(plucker_coordinates_tetrahedron[RefElement<3>::interact<1,2>(j)[i]], i);
		}
	}

};

void ComputeIntersection<Simplex<1>, Simplex<3>>::set_data(computeintersection::Simplex< 1 >& abscissa, 
                                                           computeintersection::Simplex< 3 >& tetrahedron){
	for(unsigned int i = 0; i < 4;i++){
		CI12[i].set_data(abscissa, tetrahedron[i]);
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
                    // map side (triangle) node index to tetrahedron node index
                    unsigned int node = RefElement<3>::interact<0,2>(side)[IP.idx_B()];
                    // set flag on all sides of tetrahedron connected by the node
                    for(unsigned int s=0; s < RefElement<3>::n_sides_per_node; s++)
                        CI12[RefElement<3>::interact<2,0>(node)[s]].set_computed();
                    // set topology data for object B (tetrahedron) - node
                    IP13.set_topology_B(node, IP.dim_B());
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
    
    ASSERT_LE(pocet_pruniku,2);
    
	// in the case, that line goes through vertex, but outside tetrahedron (touching vertex)
	if(pocet_pruniku == 1){
		double f_theta = IP13s[IP13s.size()-1].local_bcoords_A()[1];
        // no tolerance needed - it was already compared and normalized in 1d-2d
		if(f_theta > 1 || f_theta < 0){
			pocet_pruniku = 0;
			IP13s.pop_back();
		}

	}else if(pocet_pruniku == 2){
        
        // create shortcuts
        IntersectionPoint<1,3> 
            &a1 = IP13s[IP13s.size()-2],   // start point
            &a2 = IP13s[IP13s.size()-1];   // end point
        
        // swap the ips in the line coordinate direction (0->1 : a1->a2)        
        if(a1.local_bcoords_A()[1] > a2.local_bcoords_A()[1])
        {
//             DBGMSG("Swap.\n");
            std::swap(a1, a2);
        }
        
        // get first and second theta (coordinate of ips on the line)
        double t1, t2;
        t1 = a1.local_bcoords_A()[1];
        t2 = a2.local_bcoords_A()[1];
        
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
            arma::vec::fixed<4> interpolovane = RefElement<3>::line_barycentric_interpolation(a1.local_bcoords_B(), 
                                                                                              a2.local_bcoords_B(), 
                                                                                              a1.local_bcoords_A()[1],
                                                                                              a2.local_bcoords_A()[1], 
                                                                                              t1);
            arma::vec::fixed<2> inter({1 - t1, t1});    // barycentric coords
            a1.set_coordinates(inter,interpolovane);
//             DBGMSG("E-E 0\n");
            // set topology: node 0 of the line, tetrahedron
            a1.set_topology(0, 0, 0,3);
        }
        
        if(t1 == t2)    // if IPs are the same, then throw the second one away
        {
            pocet_pruniku = 1;
            IP13s.pop_back();
        }
        else if(t2 == 1) // interpolate IP a2
        {
            arma::vec::fixed<4> interpolovane = RefElement<3>::line_barycentric_interpolation(a1.local_bcoords_B(), 
                                                                                              a2.local_bcoords_B(), 
                                                                                              a1.local_bcoords_A()[1],
                                                                                              a2.local_bcoords_A()[1], 
                                                                                              t2);
            arma::vec::fixed<2> inter({1 - t2, t2});      // barycentric coords
            a2.set_coordinates(inter,interpolovane);
//             DBGMSG("E-E 1\n");
            // set topology: node 1 of the line, tetrahedron
            a2.set_topology(1, 0, 0,3);
        }
    }
    return pocet_pruniku;
};

void ComputeIntersection<Simplex<1>, Simplex<3>>::print_plucker_coordinates(std::ostream &os){
		os << "\tPluckerCoordinates Abscissa[0]";
		if(plucker_coordinates_abscissa_ == nullptr){
			os << "NULL" << endl;
		}else{
			os << plucker_coordinates_abscissa_;
		}

	for(unsigned int i = 0; i < 6;i++){
		os << "\tPluckerCoordinates Tetrahedron[" << i << "]";
		if(plucker_coordinates_tetrahedron[i] == nullptr){
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

// void ComputeIntersection<Simplex<1>, Simplex<3>>::set_plucker_product(double* number, unsigned int index_CI, unsigned index_edge){
// 	CI12[index_CI].set_plucker_product(number, index_edge);
// };

double* ComputeIntersection<Simplex<1>, Simplex<3>>::get_plucker_product(unsigned int side_idx, unsigned int edge_idx){
	return CI12[side_idx].get_plucker_product(edge_idx);
};






/*************************************************************************************************************
 *                                  COMPUTE INTERSECTION FOR:             2D AND 3D
 ************************************************************************************************************/
ComputeIntersection<Simplex<2>, Simplex<3>>::ComputeIntersection()
{
    plucker_coordinates_triangle_.resize(3, nullptr);
    plucker_coordinates_tetrahedron.resize(6, nullptr);
};


ComputeIntersection<Simplex<2>, Simplex<3>>::ComputeIntersection(Simplex<2> &triangle, Simplex<3> &tetrahedron)
{
    plucker_coordinates_triangle_.resize(3);
    plucker_coordinates_tetrahedron.resize(6);

    // set CI object for 1D-2D intersection 'tetrahedron edge - triangle'
	for(unsigned int i = 0; i < 6;i++){
		plucker_coordinates_tetrahedron[i] = new Plucker();
		CI12[i].set_data(tetrahedron.abscissa(i), triangle);
	}
	// set CI object for 1D-3D intersection 'triangle side - tetrahedron'
	for(unsigned int i = 0; i < 3;i++){
		plucker_coordinates_triangle_[i] = new Plucker();
		CI13[i].set_data(triangle.abscissa(i) , tetrahedron);
	}
};


void ComputeIntersection<Simplex<2>, Simplex<3>>::init(){

    // set pointers to Plucker coordinates for 1D-2D
	for(unsigned int i = 0; i < 6;i++){
		CI12[i].set_pc_abscissa(plucker_coordinates_tetrahedron[i]);
		for(unsigned int j = 0; j < 3;j++){
			CI12[i].set_pc_triangle(plucker_coordinates_triangle_[j],j);
		}
	}
	
	// set pointers to Plucker coordinates for 1D-3D
	for(unsigned int i = 0; i < 3;i++){
		CI13[i].set_pc_abscissa(plucker_coordinates_triangle_[i]);
		for(unsigned int j = 0; j < 6;j++){
			CI13[i].set_pc_tetrahedron(plucker_coordinates_tetrahedron[j],j);
		}
		CI13[i].init();
	}
};

void ComputeIntersection<Simplex<2>, Simplex<3>>::compute(computeintersection::IntersectionPolygon& local_polygon){

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
            
            local_polygon.add_ipoint(IP23);
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
				local_polygon.add_ipoint(IP23);
			}
		}
	}
};

void ComputeIntersection<Simplex<2>, Simplex<3>>::print_plucker_coordinates(std::ostream &os){
	for(unsigned int i = 0; i < 3;i++){
		os << "\tPluckerCoordinates Triangle[" << i << "]";
		if(plucker_coordinates_triangle_[i] == nullptr){
			os << "NULL" << endl;
		}else{
			os << plucker_coordinates_triangle_[i];
		}
	}
	for(unsigned int i = 0; i < 6;i++){
		os << "\tPluckerCoordinates Tetrahedron[" << i << "]";
		if(plucker_coordinates_tetrahedron[i] == nullptr){
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
