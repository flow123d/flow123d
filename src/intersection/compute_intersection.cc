/*
 *      Author: viktor
 */

#include "compute_intersection.hh"
#include "mesh/ref_element.hh"
#include "mesh/mesh.h"
#include "mesh/accessors.hh"
#include "system/system.hh"

#include "plucker.hh"
#include "intersection_point_aux.hh"
#include "intersection_aux.hh"

using namespace std;

/*************************************************************************************************************
 *                                  COMPUTE INTERSECTION FOR:             1D AND 2D
 ************************************************************************************************************/
ComputeIntersection<1,2>::ComputeIntersection()
: computed_(false)
{
    plucker_coordinates_abscissa_ = nullptr;
	plucker_coordinates_triangle_.resize(3, nullptr);
    plucker_products_.resize(3, nullptr);
}


ComputeIntersection<1,2>::ComputeIntersection(ElementAccessor<3> abscissa,
                                              ElementAccessor<3> triangle, FMT_UNUSED Mesh *mesh)
: computed_(false)
{
    ASSERT_DBG(abscissa->dim() == 1);
    ASSERT_DBG(triangle->dim() == 2);
    // in this constructor, we suppose this is the final object -> we create all data members
    plucker_coordinates_abscissa_ = new Plucker(*abscissa.node(0),
                                                *abscissa.node(1), true);
    scale_line_=plucker_coordinates_abscissa_->scale();
    
    plucker_coordinates_triangle_.resize(3);
    plucker_products_.resize(3);
    scale_triangle_=std::numeric_limits<double>::max();
    for(unsigned int side = 0; side < 3; side++){
        plucker_coordinates_triangle_[side] = new Plucker(*triangle.node(RefElement<2>::interact(Interaction<0,1>(side))[0]),
                                                          *triangle.node(RefElement<2>::interact(Interaction<0,1>(side))[1]),
                                                          true);
        scale_triangle_ = std::min( scale_triangle_, plucker_coordinates_triangle_[side]->scale());
        
        // allocate and compute new Plucker products
        plucker_products_[side] = new double((*plucker_coordinates_abscissa_)*(*plucker_coordinates_triangle_[side]));
    }
}

ComputeIntersection<1,2>::~ComputeIntersection()
{
    if(plucker_coordinates_abscissa_ != nullptr)
        delete plucker_coordinates_abscissa_;
        
    for(unsigned int side = 0; side < RefElement<2>::n_sides; side++){
        if(plucker_products_[side] != nullptr)
            delete plucker_products_[side];
        if(plucker_coordinates_triangle_[side] != nullptr)
            delete plucker_coordinates_triangle_[side];
    }
}


void ComputeIntersection<1,2>::clear_all(){
    // unset all pointers
    for(unsigned int side = 0; side < RefElement<2>::n_sides; side++)
    {
        plucker_products_[side] = nullptr;
        plucker_coordinates_triangle_[side] = nullptr;
    }
    plucker_coordinates_abscissa_ = nullptr;
}

void ComputeIntersection<1,2>::compute_plucker_products(){

    // if not already computed, compute plucker coordinates of abscissa
	plucker_coordinates_abscissa_->compute();
	scale_line_=plucker_coordinates_abscissa_->scale();

//  DBGMSG("Abscissa:\n");
//     (*abscissa_)[0].point_coordinates().print();
//     (*abscissa_)[1].point_coordinates().print();
	scale_triangle_=std::numeric_limits<double>::max();
	// if not already computed, compute plucker coordinates of triangle sides
	for(unsigned int side = 0; side < RefElement<2>::n_sides; side++){
        plucker_coordinates_triangle_[side]->compute();
		scale_triangle_ = std::min( scale_triangle_, plucker_coordinates_triangle_[side]->scale());
		
        ASSERT_DBG(plucker_products_[side]).error("Undefined plucker product.");
        if(*plucker_products_[side] == plucker_empty){
           *plucker_products_[side] = (*plucker_coordinates_abscissa_)*(*plucker_coordinates_triangle_[side]);
        }
//      DBGMSG("Plucker product = %e\n", *(plucker_products_[side]));
	}
}


int ComputeIntersection<1,2>::check_abscissa_topology(IPAux12& IP)
{
    double tol = geometry_epsilon * scale_line_;
    
    double t = RefElement<1>::bary_to_local(IP.local_bcoords_A())(0);
    
    int sign = (t < tol ? -2 : (t > 1+tol ? +2 : 0) );  // outside or inside
    if (std::abs(t) <= tol){ // left end
        sign = -1;
        IP.set_topology_A(0,0);
    }
    else if (std::abs(1-t) <= tol){ // right end
        sign = +1;
        IP.set_topology_A(1,0);
    }
    
    return sign;
}


bool ComputeIntersection<1,2>::compute_plucker(IPAux12& IP, const arma::vec3& local)
{
    // compute local barycentric coordinates of IP: see formula (3) on pg. 12 in BP VF
    // local alfa = w2/sum; local beta = w1/sum; => local barycentric coordinates in the triangle
    // where sum = w0+w1+w2

    // plucker products with correct signs according to ref_element, ordered by sides
//     double w[3] = {signed_plucker_product(0),
//                    signed_plucker_product(1),
//                    signed_plucker_product(2)};
        
    // local barycentric coordinates of IP, depends on barycentric coordinates order !!!
//     arma::vec3 local_triangle({w[2],w[1],w[0]});
//     DBGMSG("Plucker product sum = %f\n",w[0]+w[1]+w[2]);
//     local_triangle = local_triangle / (w[0]+w[1]+w[2]);
//     double w_sum = local[0] + local[1] + local[2];
//     DBGMSG("Plucker product sum = %e %e %e\n",w_sum, 1-rounding_epsilon, 1+rounding_epsilon);
    
    //assert inaccurate barycentric coordinates
    ASSERT_DBG(fabs(1.0 - local[0] - local[1] - local[2]) < geometry_epsilon)(local[0]+local[1]+local[2])
            (local[0])(local[1])(local[2]);



    arma::vec3 local_triangle({local[2],local[1],local[0]});

    // local coordinate T on the line
    // for i-th coordinate it holds: (from formula (4) on pg. 12 in BP VF)
    // T = localAbscissa= (- A(i) + ( 1 - alfa - beta ) * V0(i) + alfa * V1(i) + beta * V2 (i)) / U(i)
    // let's choose [max,i] = max {U(i)}
    arma::vec3 u = plucker_coordinates_abscissa_->get_u_vector();
    
    //find max in u in abs value:
    unsigned int i = 0; // index of maximum in u
    double max = u[0];  // maximum in u
    if(fabs((double)u[1]) > fabs(max)){ max = u[1]; i = 1;}
    if(fabs((double)u[2]) > fabs(max)){ max = u[2]; i = 2;}

    // global coordinates in triangle
    double isect_coord_i =
    local_triangle[0] * plucker_coordinates_triangle_[0]->point(RefElement<2>::normal_orientation(0))[i] +
    local_triangle[1] * plucker_coordinates_triangle_[0]->point(1-RefElement<2>::normal_orientation(0))[i] +
    local_triangle[2] * plucker_coordinates_triangle_[1]->point(RefElement<2>::normal_orientation(1))[i];

    //theta on abscissa
    double t =  (-plucker_coordinates_abscissa_->point(0)[i] + isect_coord_i)/max;
    //DebugOut() << print_var(t) << print_var(isect_coord_i) << print_var(max);
    arma::vec2 local_abscissa = {1-t, t};
    
    /*
    DBGMSG("Coordinates: line; local and global triangle\n");
    theta.print();
    local_triangle.print();
    global_triangle.print();
    */

    IP.set_topology_A(0, 1);
    IP.set_topology_B(0, 2);
    IP.set_coordinates(local_abscissa,local_triangle);
    
    return true;
}


IntersectionResult ComputeIntersection<1,2>::compute(IPAux12 &IP)
{
    compute_plucker_products();
    computed_ = true;

    //DebugOut() << "LINE: \n" << (*abscissa_)[0].point_coordinates()
    //                       << (*abscissa_)[1].point_coordinates();

    
    // convert plucker products to local coords
    arma::vec3 w = {signed_plucker_product(0),
                    signed_plucker_product(1),
                    signed_plucker_product(2)};
    double w_sum = w[0] + w[1] + w[2];
    
    unsigned int n_positive = 0;
    unsigned int n_negative = 0;
    unsigned int zero_idx_sum =0;
    //DebugOut() << print_var(std::fabs(w_sum));
    //DebugOut() << print_var(scale_line_);
    //DebugOut() << print_var(scale_triangle_);
     
    double scaled_epsilon = geometry_epsilon*scale_line_*scale_triangle_*scale_triangle_;
    if(std::fabs(w_sum) > scaled_epsilon) {
        w = w / w_sum;
        for (unsigned int i=0; i < 3; i++) {
            //DebugOut() << print_var(i) << print_var(w[i]);
            if (w[i] > geometry_epsilon) n_positive++;
            else if ( w[i] > -geometry_epsilon) zero_idx_sum+=i;
            else n_negative++;
        }

    } else {
    /* case 'w_sum == 0':
     * 1] all products are zero => n_negative=0 and n_positive=0 => it is degenerate case (coplanar case)
     * 2] at least two products are nonzero AND some of them must be negative => no intersection 
     *    (it happens when line is parallel with the triangle but not coplanar; unit test line_triangle09.msh)
     * See the IF conditions below.
     */

        for (unsigned int i=0; i < 3; i++) {
            //DebugOut().fmt("i: {} w[i]: {:g}", i, w[i]);
            if (w[i] > scaled_epsilon || w[i] < -scaled_epsilon) n_negative++;
        }
        // n_positive == 0
        //DebugOut() << print_var(n_negative);
    }

//     w.print("w");

    // any negative barycentric coordinate means, no intersection
    if (n_negative>0) return IntersectionResult::none;

    // test whether any plucker products is non-zero
    if (n_positive > 0) {
        
        compute_plucker(IP, w);
        // edge of triangle
        unsigned int non_zero_idx=0;
        if (n_positive == 2) {
            // one zero product, intersection on the zero edge
            // the zero edge index is equal to zero_idx_sum
            IP.set_topology_B(zero_idx_sum, 1);
            non_zero_idx =  (zero_idx_sum + 1) % 3;
        }
        else if (n_positive == 1) {
            // two zero products, intersection in vertex oposite to single non-zero edge
            // index of the non-zero edge is 3-zero_idx_sum
            IP.set_topology_B(RefElement<2>::oposite_node(3-zero_idx_sum), 0);
            non_zero_idx = 3-zero_idx_sum;
        }
        //DebugOut() << print_var(non_zero_idx) << print_var(signed_plucker_product(non_zero_idx));

        IntersectionResult result = signed_plucker_product(non_zero_idx) > 0 ?
                IntersectionResult::positive : IntersectionResult::negative;
        IP.set_result(result);

        return result;

    } else {
        ASSERT_DBG(IP.topology_equal(IPAux12()));   // check empty IP (none)
        IP.set_result(IntersectionResult::degenerate);  // set deg. result
        return IntersectionResult::degenerate;
    }
}

unsigned int ComputeIntersection<1,2>::compute_final(vector<IPAux12>& IP12s)
{
    IPAux12 IP;
    IntersectionResult result = compute(IP);
//     DBGVAR((int)result);
    // skip empty cases
    if(result == IntersectionResult::none) return 0;
    
    // standard case with a single intersection corner
    if(result < IntersectionResult::degenerate){
//         DBGCOUT(<< "12d plucker case\n");
        int sign = check_abscissa_topology(IP);
//         DBGVAR(sign);
        if(std::abs(sign) > 1) return 0;
        
        IP12s.push_back(IP);
        return IP12s.size();
    }
    else{
        ASSERT_DBG(result == IntersectionResult::degenerate);
//         DBGCOUT(<< "12d degenerate case, all products are zero\n");
        return compute_final_in_plane(IP12s);
    }
}


bool ComputeIntersection<1,2>::compute_degenerate(unsigned int side,
                                                  IPAux12& IP)
{
//      DBGMSG("PluckerProduct[%d]: %f\n",side, *plucker_products_[side]);
    
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
    arma::vec3 A = plucker_coordinates_abscissa_->point(0);
    // direction vector of abscissa
    arma::vec3 U = plucker_coordinates_abscissa_->get_u_vector();
    // vertex of triangle side
    arma::vec3 C = plucker_coordinates_triangle_[side]->point(RefElement<2>::normal_orientation(side));
    // direction vector of triangle side
    arma::vec3 V = plucker_coordinates_triangle_[side]->get_u_vector();
    // right hand side
    arma::vec3 K = C - A;
    // subdeterminants det_ij of the system equal minus normal vector to common plane of U and V
    // det_12 =  (-UxV)[1]
    // det_13 = -(-UxV)[2]
    // det_23 =  (-UxV)[3]
    arma::vec3 Det = -arma::cross(U,V);
    
    
//     A.print("A");
//     U.print("U");
//     C.print("C");
//     V.print("V");
//     K.print("K");
//     Det.print();
//     cout <<endl;

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
//     DBGVAR(maximum);
    //abscissa is parallel to triangle side
    //TODO What tolerance should we use here and why?
//     if(std::abs(maximum) <= std::sqrt(geometry_epsilon)) return false;
    if(std::abs(maximum) <= geometry_epsilon) return false;

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

//     DBGVAR(s);
//     DBGVAR(t);

    //TODO correct tolerance with scale; compute scale without plucker coords
    double tol = geometry_epsilon;
    // if IP is inside of triangle side
    if(t >= -tol && t <= 1+tol){

        IP.set_result(IntersectionResult::degenerate);  // set orientation as a pathologic case ( > 1)
        
        // possibly set triangle vertex {0,1,2}
        if( fabs(t) <= tol)       { IP.set_topology_B(RefElement<2>::interact(Interaction<0,1>(side))[RefElement<2>::normal_orientation(side)],0);}
        else if(fabs(1-t) <= tol) { IP.set_topology_B(RefElement<2>::interact(Interaction<0,1>(side))[1-RefElement<2>::normal_orientation(side)],0);}
        else                      { IP.set_topology_B(side,1);}   // no vertex, side side, dim = 1
        
        arma::vec2 local_abscissa({1-s, s});
        arma::vec3 local_triangle({0,0,0});

        // set local triangle barycentric coordinates according to nodes of the triangle side:
        local_triangle[RefElement<2>::interact(Interaction<0,1>(side))[RefElement<2>::normal_orientation(side)]] = 1 - t;
        local_triangle[RefElement<2>::interact(Interaction<0,1>(side))[1-RefElement<2>::normal_orientation(side)]] = t;
//      local_abscissa.print();
//      local_triangle.print();
        
        IP.set_coordinates(local_abscissa,local_triangle);
        return true; // IP found
    }
    
    return false;   // IP NOT found
}


unsigned int ComputeIntersection<1,2>::compute_final_in_plane(vector<IPAux12>& IP12s)
{
    std::vector<int> sign;
   
    for(unsigned int i = 0; i < 3;i++){
//         DBGCOUT( << "side [" << i << "]\n");
        if(IP12s.size() == 2) break;
        
        IPAux12 IP;
        if (compute_degenerate(i,IP))
        {
            uint s = check_abscissa_topology(IP);
//             cout << IP;
            // check whether we have found the same IP (e.g. vertex case)
            if(IP12s.size() > 0 && IP.topology_equal(IP12s.back())) continue;
            
            sign.push_back(s);
            IP12s.push_back(IP);
        }
    }
    if(IP12s.size() == 0) return IP12s.size();

    // in the case, that line goes through vertex, but outside tetrahedron (touching vertex)
    if(IP12s.size() == 1){
        if(std::abs(sign[0]) > 1) IP12s.pop_back(); // outside abscissa
        return IP12s.size();
    }
    
    // 2 IPs CASE:
    ASSERT_EQ_DBG(2, IP12s.size());
    
//     DBGVAR(sign[0]);
//     DebugOut() << IP12s[0];
//     DBGVAR(sign[1]);
//     DebugOut() << IP12s[1];
    
    // intersection outside of abscissa => NO intersection
    if (sign[0] == sign[1] && std::abs(sign[0]) > 1) {
        IP12s.clear();
        return 0;
    }
    
    // order IPs according to the abscissa parameter
    if(IP12s[0].local_bcoords_A()[1] > IP12s[1].local_bcoords_A()[1]){
        std::swap(IP12s[0], IP12s[1]);
        std::swap(sign[0], sign[1]);
    }
    
    // possibly cut IPs to abscissa ends and interpolate tetrahedron bcoords
    const int ip_sign[] = {-2, +2}; // states to cut
    for( unsigned int ip=0; ip<2; ip++) {
        // cut every IP to its end of the abscissa
        if (sign[ip] == ip_sign[ip]) {
            sign[ip] /=2; // -2 to -1; +2 to +1
            correct_triangle_ip_topology(double(ip), ip, IP12s);
            IP12s[ip].set_topology_A(ip, 0);
        }
    }
    
    // if IPs are the same, then throw the second one away
    if(IP12s[0].topology_equal(IP12s[1])){
//         DebugOut() << "IP same, pop up\n";
        IP12s.pop_back();
    }
    
    return IP12s.size(); 
}

void ComputeIntersection<1,2>::correct_triangle_ip_topology(
        double t, unsigned int ip, std::vector<IPAux12> &ips)
{
    arma::vec3 local_tria = RefElement<2>::line_barycentric_interpolation(
                                ips[0].local_bcoords_B(), ips[1].local_bcoords_B(),
                                ips[0].local_bcoords_A()[1], ips[1].local_bcoords_A()[1], t);
    arma::vec2 local_abscissa({1 - t, t});    // abscissa local barycentric coords
    ips[ip].set_coordinates(local_abscissa, local_tria);

    // create mask for zeros in barycentric coordinates
    // coords (*,*,*,*) -> byte bitwise xxxx
    // only least significant one byte used from the integer
    std::pair<unsigned int, unsigned int> zeros = 
        RefElement<2>::zeros_positions(ips[ip].local_bcoords_B(), geometry_epsilon);
    
    /**
     * TODO:
     * 1. Try to avoid setting topology from coords. Try to use just topology information.
     * 2. If can not be done. Use interact method to setup a map mapping 16 possible zeros positions to appropriate topology,
     *    remove topology_idx from RefElement.
     */

//     DebugOut() << "zeros: <" << zeros.first << ",  " << zeros.second << ">\n";
    switch(zeros.first)
    {
        default: ips[ip].set_topology_B(0,2);  //inside tetrahedon
                 break;
        case 1: ips[ip].set_topology_B(RefElement<2>::topology_idx<1>(zeros.second),1);
                break;
        case 2: ips[ip].set_topology_B(RefElement<2>::topology_idx<0>(zeros.second),0);
                break;
    }
}


void ComputeIntersection<1,2>::print_plucker_coordinates(std::ostream &os){
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
}




/*************************************************************************************************************
 *                                  COMPUTE INTERSECTION FOR:             2D AND 2D
 ************************************************************************************************************/
ComputeIntersection<2,2>::ComputeIntersection()
{
    plucker_coordinates_.resize(2*RefElement<2>::n_sides, nullptr);
    plucker_products_.resize(3*RefElement<2>::n_sides, nullptr);
}

ComputeIntersection<2,2>::ComputeIntersection(ElementAccessor<3> triaA,
                                              ElementAccessor<3> triaB,
                                              FMT_UNUSED Mesh *mesh)
{
    ASSERT_DBG(triaA->dim() == 2);
    ASSERT_DBG(triaB->dim() == 2);
    plucker_coordinates_.resize(2*RefElement<2>::n_sides);
    plucker_products_.resize(3*RefElement<2>::n_sides);
    
    for(unsigned int side = 0; side < RefElement<2>::n_sides; side++){
        plucker_coordinates_[side] = new Plucker(*triaA.node(RefElement<2>::interact(Interaction<0,1>(side))[0]),
                                                 *triaA.node(RefElement<2>::interact(Interaction<0,1>(side))[1]));
        plucker_coordinates_[RefElement<2>::n_sides+side]
                                   = new Plucker(*triaB.node(RefElement<2>::interact(Interaction<0,1>(side))[0]),
                                                 *triaB.node(RefElement<2>::interact(Interaction<0,1>(side))[1]));
    }

    // compute Plucker products for each pair triangle A side and triangle B side
    for(unsigned int p = 0; p < 3*RefElement<2>::n_sides; p++){
        plucker_products_[p] = new double(plucker_empty);
    }
}

ComputeIntersection<2,2>::~ComputeIntersection()
{
    // unset pointers:
    for(unsigned int side = 0; side <  2*RefElement<2>::n_sides; side++)
        CI12[side].clear_all();
    
    // then delete objects:
    for(unsigned int side = 0; side < 2*RefElement<2>::n_sides; side++){
        if(plucker_coordinates_[side] != nullptr)
            delete plucker_coordinates_[side];
    }
    
    for(unsigned int p = 0; p < 3*RefElement<2>::n_sides; p++){
        if(plucker_products_[p] != nullptr)
            delete plucker_products_[p];
    }
}

void ComputeIntersection<2,2>::clear_all()
{
    // unset all pointers
    for(unsigned int side = 0; side < 2*RefElement<2>::n_sides; side++)
    {
        plucker_coordinates_[side] = nullptr;
    }
    for(unsigned int p = 0; p < 3*RefElement<2>::n_sides; p++){
        plucker_products_[p] = nullptr;
    }
}

void ComputeIntersection<2,2>::init(){

//     DBGMSG("init\n");
    for(unsigned int i = 0; i <  RefElement<2>::n_sides; i++){
        // set side A vs triangle B
        for(unsigned int j = 0; j <  RefElement<2>::n_sides; j++)
            CI12[i].set_pc_triangle(plucker_coordinates_[3+j], j);  // set triangle B
        // set side of triangle A
        CI12[i].set_pc_abscissa(plucker_coordinates_[i]);
        
        // set side B vs triangle A
        for(unsigned int j = 0; j <  RefElement<2>::n_sides; j++)
            CI12[RefElement<2>::n_sides + i].set_pc_triangle(plucker_coordinates_[j], j); // set triangle A
        // set side of triangle B
        CI12[RefElement<2>::n_sides + i].set_pc_abscissa(plucker_coordinates_[RefElement<2>::n_sides + i]);
        
        // set plucker products
        for(unsigned int j = 0; j <  RefElement<2>::n_sides; j++)
        {
            //for A[i]_B set pp. A[i] x B[j]
            CI12[i].set_plucker_product(get_plucker_product(i,j),j);
            //for B[i]_A set pp. A[j] x B[i]
            CI12[RefElement<2>::n_sides + i].set_plucker_product(get_plucker_product(j,i),j);
        }
        
    }
}


unsigned int ComputeIntersection<2,2>::compute(IntersectionAux<2,2>& intersection)
{
    // final intersection points
    std::vector<IPAux22> &IP22s = intersection.points();
    // temporary vector for lower dimensional IPs
    std::vector<IPAux12> IP12s;
    IP12s.reserve(2);
    unsigned int ip_coutner = 0;

    // loop over CIs (side vs triangle): [A0_B, A1_B, A2_B, B0_A, B1_A, B2_A].
    for(unsigned int i = 0; i < 2*RefElement<2>::n_sides && ip_coutner < 2; i++){
        if(!CI12[i].is_computed()) // if not computed yet
        {
//             DBGVAR(i);
            if(CI12[i].compute_final(IP12s) == 0) continue;
            
            unsigned int triangle_side = i%3; //i goes from 0 to 5 -> i%3 = 0,1,2,0,1,2
            
            for(IPAux12 &IP : IP12s)
            {
                IPAux22 IP22(IP.switch_objects(), triangle_side);   // swicth dim 12 -> 21; interpolate dim 21 -> 22
                
                if(i < 3){
                    //switch back to keep order of triangles [A,B]
                    IP22 = IP22.switch_objects();
            
                    if( IP22.dim_A() == 0 ) // if IP is vertex of triangle
                    {
//                         DBGCOUT("set_node A\n");
                        // we are on line of the triangle A, and IP.idx_A contains local node of the line
                        // we know vertex index
                        
                        // set flag on all sides of tetrahedron connected by the node
                        for(unsigned int s=0; s < RefElement<2>::n_sides_per_node; s++)
                            CI12[RefElement<2>::interact(Interaction<1,0>(IP22.idx_A()))[s]].set_computed();
                    }
                    if( IP22.dim_B() == 0 ) // if IP is vertex of triangle
                    {
//                         DBGCOUT("set_node B\n");
                        // set flag on both sides of triangle connected by the node
                        for(unsigned int s=0; s < RefElement<2>::n_sides_per_node; s++)
                            CI12[3 + RefElement<2>::interact(Interaction<1,0>(IP22.idx_B()))[s]].set_computed();
                    }
                    else if( IP22.dim_B() == 1 ) // if IP is vertex of triangle
                    {
//                         DBGCOUT("set line B\n");
                        // set flag on both sides of triangle connected by the node
                        CI12[3 + IP22.idx_B()].set_computed();
                    }
                }
                else if( IP22.dim_B() == 0 ) // if IP is vertex of triangle B (triangles switched!!!  A <-> B)
                {
                    //we do not need to look back to triangle A, because if IP was at vertex, we would know already
//                     DBGCOUT("set_node B\n");
                    // we are on line of the triangle, and IP.idx_A contains local node of the line
                    // we know vertex index
                    
                    // set flag on both sides of triangle connected by the node
                    for(unsigned int s=0; s < RefElement<2>::n_sides_per_node; s++)
                        CI12[3 + RefElement<2>::interact(Interaction<1,0>(IP22.idx_B()))[s]].set_computed();
                }
//                 DBGCOUT( << IP22);
                ip_coutner++;
                IP22s.push_back(IP22);
            }
            IP12s.clear();
        }
    }

    return ip_coutner;
}

// void ComputeIntersection< 2, 2>::correct_triangle_ip_topology(IntersectionPointAux<2,2>& ip)
// {
//     // create mask for zeros in barycentric coordinates
//     // coords (*,*,*,*) -> byte bitwise xxxx
//     // only least significant one byte used from the integer
//     unsigned int zeros = 0;
//     unsigned int n_zeros = 0;
//     for(char i=0; i < 3; i++){
//         if(std::fabs(ip.local_bcoords_B()[i]) < geometry_epsilon)
//         {
//             zeros = zeros | (1 << (2-i));
//             n_zeros++;
//         }
//     }
//     
//     switch(n_zeros)
//     {
//         default: ip.set_topology_B(0,2);  //inside triangle
//                  break;
//         case 1: ip.set_topology_B(RefElement<2>::topology_idx<1>(zeros),2);
//                 break;
//         case 2: ip.set_topology_B(RefElement<2>::topology_idx<0>(zeros),1);
//                 break;
//     }
// }


void ComputeIntersection<2,2>::print_plucker_coordinates(std::ostream &os){

    for(unsigned int i = 0; i < RefElement<2>::n_lines; i++){
        os << "\tPluckerCoordinates Triangle A[" << i << "]";
        if(plucker_coordinates_[i] == nullptr){
            os << "NULL" << endl;
        }else{
            os << plucker_coordinates_[i];
        }
    }
    for(unsigned int i = 0; i < RefElement<2>::n_lines; i++){
        os << "\tPluckerCoordinates Triangle B[" << i << "]";
        if(plucker_coordinates_[RefElement<2>::n_lines + i] == nullptr){
            os << "NULL" << endl;
        }else{
            os << plucker_coordinates_[RefElement<2>::n_lines + i];
        }
    }
}

void ComputeIntersection<2,2>::print_plucker_coordinates_tree(std::ostream &os){
    os << "ComputeIntersection<2,2> Plucker Coordinates Tree:" << endl;
        print_plucker_coordinates(os);
        for(unsigned int i = 0; i < 6;i++){
            os << "ComputeIntersection<1,2>["<< i <<"] Plucker Coordinates:" << endl;
            CI12[i].print_plucker_coordinates(os);
        }
}


/*************************************************************************************************************
 *                                  COMPUTE INTERSECTION FOR:             1D AND 3D
 ************************************************************************************************************/
ComputeIntersection<1,3>::ComputeIntersection()
{
    plucker_coordinates_abscissa_ = nullptr;
    plucker_coordinates_tetrahedron.resize(6, nullptr);
    plucker_products_.resize(6, nullptr);
}

ComputeIntersection<1,3>::ComputeIntersection(ElementAccessor<3> abscissa,
                                              ElementAccessor<3> tetrahedron,
                                              FMT_UNUSED Mesh *mesh)
{
    ASSERT_DBG(abscissa->dim() == 1);
    ASSERT_DBG(tetrahedron->dim() == 3);
    
    plucker_coordinates_abscissa_ = new Plucker(*abscissa.node(0), *abscissa.node(1));
    plucker_coordinates_tetrahedron.resize(6);
    plucker_products_.resize(6);
    
    for(unsigned int line = 0; line < RefElement<3>::n_lines; line++){
        plucker_coordinates_tetrahedron[line] = new Plucker(*tetrahedron.node(RefElement<3>::interact(Interaction<0,1>(line))[0]),
                                                            *tetrahedron.node(RefElement<3>::interact(Interaction<0,1>(line))[1]));
        // compute Plucker products (abscissa X tetrahedron line)
        plucker_products_[line] = new double(plucker_empty);
    }
}

ComputeIntersection<1,3>::~ComputeIntersection()
{
    // unset pointers:
    for(unsigned int side = 0; side <  RefElement<3>::n_sides; side++)
        CI12[side].clear_all();
    
    // then delete objects:
    if(plucker_coordinates_abscissa_ != nullptr)
        delete plucker_coordinates_abscissa_;
    
    for(unsigned int line = 0; line < RefElement<3>::n_lines; line++){
        if(plucker_products_[line] != nullptr)
            delete plucker_products_[line];
        if(plucker_coordinates_tetrahedron[line] != nullptr)
            delete plucker_coordinates_tetrahedron[line];
    }
}

void ComputeIntersection<1,3>::clear_all()
{
    // unset all pointers
    for(unsigned int side = 0; side < RefElement<3>::n_lines; side++)
    {
        plucker_products_[side] = nullptr;
        plucker_coordinates_tetrahedron[side] = nullptr;
    }
    plucker_coordinates_abscissa_ = nullptr;
}

void ComputeIntersection<1,3>::init(){

	for(unsigned int side = 0; side <  RefElement<3>::n_sides; side++){
		for(unsigned int line = 0; line < RefElement<3>::n_lines_per_side; line++){
			CI12[side].set_pc_triangle(plucker_coordinates_tetrahedron[
                                            RefElement<3>::interact(Interaction<1,2>(side))[line]], line);
            
            CI12[side].set_plucker_product(
                plucker_products_[RefElement<3>::interact(Interaction<1,2>(side))[line]],
                line);
		}
		CI12[side].set_pc_abscissa(plucker_coordinates_abscissa_);
	}  
}


unsigned int ComputeIntersection<1,3>::compute(IntersectionAux< 1, 3 >& intersection)
{
    return compute(intersection.i_points_);
}

unsigned int ComputeIntersection<1,3>::compute(std::vector<IPAux> &IP13s){
	ASSERT_EQ_DBG(0, IP13s.size());
    std::vector<int> sign;

   // loop over faces of tetrahedron
	for(unsigned int face = 0; face < RefElement<3>::n_sides && IP13s.size() < 2; face++){
        
		if  (CI12[face].is_computed()) continue;
        
        IntersectionPointAux<1,2> IP12;
	    IntersectionResult result = CI12[face].compute(IP12);
        //DebugOut() << print_var(face) << print_var(int(result)) << "1d-3d";

		if (int(result) < int(IntersectionResult::degenerate) ) {
            // fix topology of A (abscissa) in IP12 and save sign
            sign.push_back(CI12[face].check_abscissa_topology(IP12));
            
            IPAux IP13(IP12, face);
            //DebugOut() << IP13;
            IP13s.push_back(IP13);
            
            // set the 'computed' flag on the connected sides by IP
            if(IP13.dim_B() == 0) // IP is vertex of triangle
            {
                // set flag on all sides of tetrahedron connected by the node
                for(unsigned int node_face : RefElement<3>::interact(Interaction<2,0>(IP13.idx_B())))
                    CI12[node_face].set_computed();
            }
            else if(IP13.dim_B() == 1) // IP is on edge of triangle
            {
                for(unsigned int edge_face : RefElement<3>::interact(Interaction<2,1>(IP13.idx_B())))
                    CI12[edge_face].set_computed();
            }
		}
	}
	if (IP13s.size() == 0) return IP13s.size();

	// in the case, that line goes through vertex, but outside tetrahedron (touching vertex)
	if(IP13s.size() == 1){
        if(std::abs(sign[0]) > 1) IP13s.pop_back(); // outside abscissa
        return IP13s.size();
	}
    
    // 2 IPs CASE:
    ASSERT_EQ_DBG(2, IP13s.size());
    
    // intersection outside of abscissa => NO intersection
    if (sign[0] == sign[1] && std::abs(sign[0]) > 1) {
        IP13s.clear();
        return 0;
    }
    
    // order IPs according to the abscissa parameter
    if(IP13s[0].local_bcoords_A()[1] > IP13s[1].local_bcoords_A()[1]){
        std::swap(IP13s[0], IP13s[1]);
        std::swap(sign[0], sign[1]);
    }
    
    // possibly cut IPs to abscissa ends and interpolate tetrahedron bcoords
    const int ip_sign[] = {-2, +2}; // states to cut
    for( unsigned int ip=0; ip<2; ip++) {
        // cut every IP to its end of the abscissa
        if (sign[ip] == ip_sign[ip]) {
            sign[ip] /=2; // -2 to -1; +2 to +1
            correct_tetrahedron_ip_topology(double(ip), ip, IP13s);
            IP13s[ip].set_topology_A(ip, 0);
        }
    }

    // if IPs are the same, then throw the second one away
    if(IP13s[0].topology_equal(IP13s[1])){
        IP13s.pop_back();
    }
    
    return IP13s.size();
}

void ComputeIntersection<1,3>::correct_tetrahedron_ip_topology(
        double t, unsigned int ip, std::vector<IPAux> &ips)
{
    arma::vec4 local_tetra = RefElement<3>::line_barycentric_interpolation(
                                ips[0].local_bcoords_B(), ips[1].local_bcoords_B(),
                                ips[0].local_bcoords_A()[1], ips[1].local_bcoords_A()[1], t);
    arma::vec2 local_abscissa({1 - t, t});    // abscissa local barycentric coords
    ips[ip].set_coordinates(local_abscissa, local_tetra);

    // create mask for zeros in barycentric coordinates
    // coords (*,*,*,*) -> byte bitwise xxxx
    // only least significant one byte used from the integer
    unsigned int zeros = 0;
    unsigned int n_zeros = 0;
    for(char i=0; i < 4; i++){
        if(std::fabs(ips[ip].local_bcoords_B()[i]) < geometry_epsilon)
        {
            zeros = zeros | (1 << i);
            n_zeros++;
        }
    }
    
    /**
     * TODO:
     * 1. Try to avoid setting topology from coords. Try to use just topology information.
     * 2. If can not be done. Use interact method to setup a map mapping 16 possible zeros positions to appropriate topology,
     *    remove topology_idx from RefElement.
     */

    switch(n_zeros)
    {
        default: ips[ip].set_topology_B(0,3);  //inside tetrahedon
                 break;
        case 1: ips[ip].set_topology_B(RefElement<3>::topology_idx<2>(zeros),2);
                break;
        case 2: ips[ip].set_topology_B(RefElement<3>::topology_idx<1>(zeros),1);
                break;
        case 3: ips[ip].set_topology_B(RefElement<3>::topology_idx<0>(zeros),0);
                break;
    }
}


void ComputeIntersection<1,3>::print_plucker_coordinates(std::ostream &os){
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
}

void ComputeIntersection<1,3>::print_plucker_coordinates_tree(std::ostream &os){
	os << "ComputeIntersection<1,3> Plucker Coordinates Tree:" << endl;
		print_plucker_coordinates(os);
		for(unsigned int i = 0; i < 4;i++){
			os << "ComputeIntersection<1,2>["<< i <<"] Plucker Coordinates:" << endl;
			CI12[i].print_plucker_coordinates(os);
		}
}



/*************************************************************************************************************
 *                                  COMPUTE INTERSECTION FOR:             2D AND 3D
 ************************************************************************************************************/
ComputeIntersection<2,3>::ComputeIntersection()
: no_idx(100),
s3_dim_starts({0, 4, 10, 14}), // vertices, edges, faces, volume
s2_dim_starts({15, 18, 21}),   // vertices, sides, surface
object_next(22, no_idx),       // 4 vertices, 6 edges, 4 faces, 1 volume, 3 corners, 3 sides, 1 surface; total 22
on_faces(_on_faces())
 {

    plucker_coordinates_triangle_.resize(3, nullptr);
    plucker_coordinates_tetrahedron.resize(6, nullptr);
    plucker_products_.resize(3*6, nullptr);
}


ComputeIntersection<2,3>::ComputeIntersection(ElementAccessor<3> triangle,
                                              ElementAccessor<3> tetrahedron,
                                              Mesh *mesh)
: ComputeIntersection()
{
    mesh_ = mesh;
    plucker_coordinates_triangle_.resize(3);
    plucker_coordinates_tetrahedron.resize(6);

    // set CI object for 1D-2D intersection 'tetrahedron edge - triangle'
	for(unsigned int i = 0; i < RefElement<3>::n_lines; i++){
		plucker_coordinates_tetrahedron[i] = new Plucker(*tetrahedron.node(RefElement<3>::interact(Interaction<0,1>(i))[0]),
                                                         *tetrahedron.node(RefElement<3>::interact(Interaction<0,1>(i))[1]));
	}
	// set CI object for 1D-3D intersection 'triangle side - tetrahedron'
	for(unsigned int i = 0; i < RefElement<2>::n_lines;i++){
		plucker_coordinates_triangle_[i] = new Plucker(*triangle.node(RefElement<2>::interact(Interaction<0,1>(i))[0]),
                                                       *triangle.node(RefElement<2>::interact(Interaction<0,1>(i))[1]));
	}
	
	// compute Plucker products (triangle side X tetrahedron line)
	// order: triangle sides X tetrahedron lines:
	// TS[0] X TL[0..6]; TS[1] X TL[0..6]; TS[1] X TL[0..6]
	unsigned int np = RefElement<2>::n_sides *  RefElement<3>::n_lines;
	plucker_products_.resize(np, nullptr);
    for(unsigned int line = 0; line < np; line++){
        plucker_products_[line] = new double(plucker_empty);
        
    }
}

ComputeIntersection<2,3>::~ComputeIntersection()
{
    // unset pointers:
    for(unsigned int triangle_side = 0; triangle_side < RefElement<2>::n_sides; triangle_side++)
        CI13[triangle_side].clear_all();
        
    for(unsigned int line = 0; line < RefElement<3>::n_lines; line++)
        CI12[line].clear_all();
    
    // then delete objects:
    unsigned int np = RefElement<2>::n_sides *  RefElement<3>::n_lines;
    for(unsigned int line = 0; line < np; line++){
            if(plucker_products_[line] != nullptr)
                delete plucker_products_[line];
    }
    
    for(unsigned int i = 0; i < RefElement<3>::n_lines;i++){
        if(plucker_coordinates_tetrahedron[i] != nullptr)
            delete plucker_coordinates_tetrahedron[i];
    }
    for(unsigned int i = 0; i < RefElement<2>::n_sides;i++){
        if(plucker_coordinates_triangle_[i] != nullptr)
            delete plucker_coordinates_triangle_[i];
    }
}

void ComputeIntersection<2,3>::init(){

    // set pointers to Plucker coordinates for 1D-2D
    // set pointers to Plucker coordinates for 1D-3D
    // distribute Plucker products - CI13
    for(unsigned int triangle_side = 0; triangle_side < RefElement<2>::n_sides; triangle_side++){
        for(unsigned int line = 0; line < RefElement<3>::n_lines; line++){
            CI13[triangle_side].set_plucker_product(
                    plucker_products_[triangle_side * RefElement<3>::n_lines + line],
                    line);
            CI12[line].set_plucker_product(
                    plucker_products_[triangle_side * RefElement<3>::n_lines + line],
                    triangle_side);
            
            CI13[triangle_side].set_pc_tetrahedron(plucker_coordinates_tetrahedron[line],line);
            CI12[line].set_pc_triangle(plucker_coordinates_triangle_[triangle_side],triangle_side);
        }
        CI13[triangle_side].set_pc_abscissa(plucker_coordinates_triangle_[triangle_side]);
        CI13[triangle_side].init();
    }
    
    // set pointers to Plucker coordinates for 1D-2D
    for(unsigned int line = 0; line < RefElement<3>::n_lines; line++)
        CI12[line].set_pc_abscissa(plucker_coordinates_tetrahedron[line]);

}



bool ComputeIntersection<2,3>::have_backlink(uint i_obj) {
    ASSERT_LT_DBG(i_obj, object_next.size());
    unsigned int ip = object_next[i_obj];
    if (ip == no_idx) return false;
    ASSERT_LT_DBG(ip, IP_next.size());
    return IP_next[ip] == i_obj;
}

/**
 * Set links: obj_before -> IP -> obj_after
 * if obj_after have null successor, set obj_after -> IP (backlink)
 */
void ComputeIntersection<2,3>::set_links(uint obj_before_ip, uint ip_idx, uint obj_after_ip) {
    if (have_backlink(obj_after_ip)) {
        // target object is already target of other IP, so it must be source object
        std::swap(obj_before_ip, obj_after_ip);
    }
    //DebugOut().fmt("before: {} ip: {} after: {}\n", obj_before_ip, ip_idx, obj_after_ip );
    ASSERT_DBG( ! have_backlink(obj_after_ip) )
        (mesh_->element_accessor(intersection_->component_ele_idx()).idx())
        (mesh_->element_accessor(intersection_->bulk_ele_idx()).idx())
        (obj_before_ip)(ip_idx)(obj_after_ip); // at least one could be target object
    object_next[obj_before_ip] = ip_idx;
    IP_next.push_back( obj_after_ip);
    if (object_next[obj_after_ip] == no_idx) {
        object_next[obj_after_ip] = ip_idx;
    }
}



void ComputeIntersection<2,3>::compute(IntersectionAux< 2 , 3  >& intersection)
{
    intersection_= &intersection;
    //DebugOut().fmt("2d ele: {} 3d ele: {}\n",
    //        intersection.component_ele_idx(),
    //        intersection.bulk_ele_idx());

    IP23_list.clear();
    IP_next.clear();
    std::fill(object_next.begin(), object_next.end(), no_idx);
	std::vector<IPAux13> IP13s;

	std::array<bool, 6> edge_touch={false,false, false, false, false, false};
    unsigned int object_before_ip, object_after_ip;

    //unsigned int last_triangle_vertex=30; // no vertex at last IP

	// pass through the ccwise oriented sides in ccwise oriented order
	// How to make this in independent way?
	// Move this into RefElement?
	std::vector<unsigned int> side_cycle_orientation = { 0, 0, 1};
	std::vector<unsigned int> cycle_sides = {0, 2, 1};

	// TODO:
	// better mechanism for detecting vertex duplicities, do no depend on cyclic order of sides
	// still need cyclic orientation
	for(unsigned int _i_side = 0; _i_side < RefElement<2>::n_lines; _i_side++) {    // go through triangle lines
	    unsigned int i_side = cycle_sides[_i_side];
	    IP13s.clear();
        CI13[ i_side ].compute(IP13s);
        ASSERT_DBG(IP13s.size() < 3);
        if (IP13s.size() == 0) continue;
        for(unsigned int _ip=0; _ip < IP13s.size(); _ip++) {
            //int ip_topo_position = _ip*2-1; // -1 (incoming ip), +1 (outcoming ip), 0 both

            // fix order of IPs
            unsigned int ip = (side_cycle_orientation[_i_side] + _ip) % IP13s.size();

//             DebugOut().fmt("rside: {} cside: {} rip: {} cip: {}", _i_side, i_side, _ip, ip);

            // convert from 13 to 23 IP
            IntersectionPointAux<3,1> IP31 = IP13s[ip].switch_objects();   // switch idx_A and idx_B and coords
            IntersectionPointAux<3,2> IP32(IP31, i_side);    // interpolation uses local_bcoords_B and given idx_B
            IPAux23 IP23 = IP32.switch_objects(); // switch idx_A and idx_B and coords back
            //DebugOut() << IP;

            // Tracking info
            unsigned int tetra_object = s3_dim_starts[IP23.dim_B()] + IP23.idx_B();
            unsigned int side_object = s2_dim_starts[1] + i_side;

            object_before_ip = tetra_object;
            object_after_ip = side_object;


            // IP is vertex of triangle,
            if( IP23.dim_A() == 0 && IP23.dim_B() == 3)
            {
                // we are on line of the triangle, and IP.idx_A contains local node of the line
                // E-E, we know vertex index
                object_before_ip = s2_dim_starts[0]+IP23.idx_A();
            }// else current_triangle_vertex=3+IP23_list.size(); // no vertex, and unique

            // side of triangle touching  S3, in vertex or in edge
            if (IP13s.size() == 1 ) {
                if (IP23.dim_B() == 0) {
                    continue; // skip, S3 vertices are better detected in phase 2
                }
                if (IP23.dim_A() == 0) { // vertex of triangle
                    object_before_ip = tetra_object;
                    object_after_ip = s2_dim_starts[0]+IP23.idx_A();

                    // source vertex of the side vector (oriented CCwise)
                    //if ( (IP.idx_A()+side_cycle_orientation[_i_side])%2 == 0)
                    //    std::swap(object_before_ip, object_after_ip);
                } else {
                    // touch in edge

                    //continue;
                    ASSERT_EQ_DBG(IP23.dim_B(), 1);
                    edge_touch[IP23.idx_B()]=true;
                    std::swap(object_before_ip, object_after_ip);
                }
            }

            IP23_list.push_back(IP23);

            unsigned int ip_idx = IP23_list.size()-1;
            //DebugOut().fmt("before: {} after: {} ip: {}\n", object_before_ip, object_after_ip, ip_idx);
            ASSERT_EQ_DBG(IP23_list.size(),  IP_next.size()+1);
            set_links(object_before_ip, ip_idx, object_after_ip);
        }

    }

    // now we have at most single true degenerate IP in IP23
    // TODO:
    // - deal with degenerate IPs in the edge-trinagle phase
    // - remove degenerate_list (just count degen points)
    // - remove check for duplicities in final list copy
    // - add more comment to individual cases in order to be sure that any case in particular branch is
    //   treated right


    //TODO: cannot the two cycles be merged in one now?
    // and instead of processed_edge use CI12[tetra_edge].set_computed();
    // and I think we don't have to generate dummy IP12s
    // so we don't have to have IP12s vector here
    IP12s_.clear();
    // S3 Edge - S2 intersections; collect all signs, make dummy intersections
	for(unsigned int tetra_edge = 0; tetra_edge < 6; tetra_edge++) {
        IPAux12 IP12;
        CI12[tetra_edge].compute(IP12);
	    // IntersectionResult result = CI12[tetra_edge].compute(IP12);
	    // DebugOut() << print_var(tetra_edge) << print_var(int(result));
        // in degenerate case: IP12 is empty with degenerate result
        IP12s_.push_back(IP12);
	}
	vector<uint> processed_edge(6, 0);
	FacePair face_pair;
	for(unsigned int tetra_edge = 0; tetra_edge < 6;tetra_edge++) {
	    if (! processed_edge[tetra_edge]) {
//             DBGVAR(tetra_edge);
	        IPAux12 &IP12 = IP12s_[tetra_edge];
            if(IP12.result() >= IntersectionResult::degenerate) continue;
            
            int sign = CI12[tetra_edge].check_abscissa_topology(IP12);
//             DBGVAR(sign);
            if(std::abs(sign) > 1) continue;
            
            IPAux23 IP23(IP12s_[tetra_edge].switch_objects(), tetra_edge);
            
            const uint edge_dim = IP23.dim_B();
            const uint i_edge = IP23.idx_B();
	        ASSERT_LT_DBG(edge_dim, 2);

            if ( edge_dim == 0) {
                face_pair = vertex_faces(i_edge);
                // mark edges coincident with the vertex
                for( uint ie : RefElement<3>::interact(Interaction<1,0>(i_edge)) )
                    processed_edge[ie] = 1;
            }
            else
                face_pair = edge_faces(i_edge);
                
            //DebugOut() << print_var(face_pair[0])<< print_var(face_pair[1]);

            IP23_list.push_back(IP23);
            unsigned int ip_idx = IP23_list.size()-1;

            unsigned int s3_object = s3_dim_starts[edge_dim] + i_edge;

            //DebugOut() << print_var(edge_touch[i_edge]) << print_var(have_backlink(s3_object));
	        if (IP23.dim_A() < 2
	                && (! edge_touch[i_edge])
	                && object_next[s3_object] != no_idx) { // boundary of S2, these ICs are duplicate

	            if ( have_backlink(s3_object) ) {
	                set_links(s3_object, ip_idx, face_pair[1]);
	            } else {
	                set_links(face_pair[0], ip_idx, s3_object);
	            }

	        } else { // interior of S2, just use the face pair
	                //DebugOut() << print_var(face_pair[0])<< print_var(face_pair[1]);
	                set_links(face_pair[0], ip_idx, face_pair[1]);

	                if ( have_backlink(s3_object) ) {
	                    object_next[s3_object]=ip_idx;
	                }
	        }
	    }
	}


    // Return IPs in correct order and remove duplicates
	ASSERT_EQ(0, intersection.size());

    if (IP23_list.size() == 0) return; // empty intersection

    // detect first IP, this needs to be done only in the case of
    // point or line intersections, where IPs links do not form closed cycle
    // Possibly we do this only if we detect such case through previous phases.
    vector<char> have_predecessor(IP23_list.size(), 0);
    for(auto obj : IP_next) {
        ASSERT_LT_DBG(obj, object_next.size());
        unsigned int ip = object_next[obj];
        if (ip < IP_next.size()) have_predecessor[ip]=1;
    }
    unsigned int ip_init=0;
    for(unsigned int i=0; i< IP23_list.size(); i++) if (! have_predecessor[i]) ip_init=i;

    arma::uvec::fixed<RefElement<3>::n_sides> ips_face_counter; 
    ips_face_counter.zeros();
    
    // regular case, merge duplicit IPs
    unsigned int ip=ip_init;
    ASSERT_EQ_DBG(IP_next.size(), IP23_list.size());
    intersection.points().push_back(IP23_list[ip]);
    //DebugOut() << print_var(ip) << IP23_list[ip];
    while (1)  {
        //DebugOut() << print_var(ip) << IP23_list[ip];

        unsigned int object = IP_next[ip];
        //IP_next[ip]=no_idx;
        ASSERT_LT_DBG(object, object_next.size());
        ip = object_next[object];
        object_next[object]=no_idx;
        if ((ip == no_idx)) break;
        ASSERT_LT_DBG(ip, IP_next.size());

        if ( ! IP23_list[ip].topology_equal(intersection.points().back()) ) {
            IPAux23 &IP = IP23_list[ip];
            //DebugOut() << print_var(ip) << IP23_list[ip];
            intersection.points().push_back(IP);
            if(IP.dim_B() < 3)
                ips_face_counter += on_faces[IP.dim_B()][IP.idx_B()];
        }
    }
    
    if (intersection.points().size() == 1) return;

    if (IP23_list[ip_init].topology_equal(intersection.points().back()) )
        intersection.points().pop_back();
    
    if (intersection.points().size() > 2){
        for(uint i=0; i< ips_face_counter.n_elem; i++){
            if(ips_face_counter(i) >= intersection.points().size()){ // all IPs lie in a single face
//                 DBGCOUT(<< "all IPs in face: " << i <<"\n");
//                 ips_face_counter.print("face_counter");
                intersection.set_ips_in_face(i);
            }
        }
    }
}

auto ComputeIntersection<2,3>::edge_faces(uint i_edge)-> FacePair
{
    auto &line_faces=RefElement<3>::interact(Interaction<2,1>(i_edge));
    unsigned int ip_ori = (unsigned int)(IP12s_[i_edge].result());
    ASSERT_DBG(ip_ori < 2); // no degenerate case

    // RefElement returns edge faces in clockwise order (edge pointing to us)
    // negative ip sign (ori 0) = faces counter-clockwise
    // positive ip sign (ori 1) = faces clockwise
    return { s3_dim_starts[2] + line_faces[1-ip_ori], s3_dim_starts[2] + line_faces[ip_ori] };
}

auto ComputeIntersection<2,3>::vertex_faces(uint i_vertex)-> FacePair
{
    // vertex edges clockwise
    const IdxVector<3> &vtx_edges = RefElement<3>::interact(Interaction<1,0>(i_vertex));
    std::array<unsigned int, 3> n_ori, sum_idx;
    n_ori.fill(0);
    sum_idx.fill(0);
    for(unsigned int ie=0; ie <3; ie++) {
        unsigned int edge_ip_ori = (unsigned int)(IP12s_[ vtx_edges[ie]].result());
        if (RefElement<3>::interact(Interaction<0,1>(vtx_edges[ie]))[0] != i_vertex
            && edge_ip_ori!= int(IntersectionResult::degenerate) )
            edge_ip_ori = (edge_ip_ori +1)%2;
        //ASSERT_LT_DBG(edge_ip_ori, 3)(ie);
        if (edge_ip_ori == 3) edge_ip_ori=2; // none as degenerate
        n_ori[edge_ip_ori]++;
        sum_idx[edge_ip_ori]+=ie;
    }
    unsigned int n_degen = n_ori[ int(IntersectionResult::degenerate) ];
    unsigned int sum_degen = sum_idx[ int(IntersectionResult::degenerate) ];
    unsigned int n_positive = n_ori[ int(IntersectionResult::positive) ];
    unsigned int n_negative= n_ori[ int(IntersectionResult::negative) ];
    //DebugOut().fmt("nd: {} sd: {} np: {} nn: {}", n_degen, sum_degen, n_positive, n_negative);
    if ( n_degen == 2 ) {
        // S2 plane match a face of S3, we treat degenerated edges as the faces
        // incident with the single regualr edge.

        unsigned int i_edge = 3 - sum_degen; // regular edge index
        FacePair pair = edge_faces(vtx_edges[i_edge]);
        auto &vtx_faces = RefElement<3>::interact(Interaction<2,0>(i_vertex));
        // replace faces by edges
        if (pair[0] == s3_dim_starts[2] + vtx_faces[(i_edge+1)%3])
            return { s3_dim_starts[1] + (i_edge+2)%3,  s3_dim_starts[1] + (i_edge+1)%3 };
        else
            return { s3_dim_starts[1] + (i_edge+1)%3,  s3_dim_starts[1] + (i_edge+2)%3 };

    } else if (n_degen == 1) {
        // One edge in S2 plane.
        unsigned int i_edge = sum_degen;
        ASSERT( n_positive + n_negative == 2);
        if ( n_positive == 1) {
            // opposite signs, S2 plane cuts S3

            FacePair pair = edge_faces(vtx_edges[(i_edge+1)%3]);
            unsigned int face = RefElement<3>::interact(Interaction<2,0>(i_vertex))[i_edge];
            // assign edges to faces
            //DebugOut().fmt("vtx: {} edg: {} face: {}", i_vertex, i_edge, face);
            if (pair[0] == s3_dim_starts[2] + face)
                return { s3_dim_starts[2] + face, s3_dim_starts[1] + vtx_edges[i_edge]};
            else
                return { s3_dim_starts[1] + vtx_edges[i_edge], s3_dim_starts[2] + face };
        } else {
            // same signs; S2 plane touch S3 vertex and a single edge
            //DebugOut() << "Touch in edge.";
            // same signs S2 plane touchs S3
            ASSERT(n_positive == 0 || n_positive== 2);
            return { s3_dim_starts[0]+i_vertex, s3_dim_starts[1] + vtx_edges[i_edge]};
        }


    } else {
        ASSERT(n_degen == 0);
        ASSERT( n_positive + n_negative == 3);

        if (n_positive == 1) {
            unsigned int i_edge = sum_idx[ int(IntersectionResult::positive) ];
            return edge_faces(vtx_edges[i_edge]);
        } else if (n_negative == 1) {
            unsigned int i_edge = sum_idx[ int(IntersectionResult::negative) ];
            return edge_faces(vtx_edges[i_edge]);
        } else {
            // S2 touch vertex of S3 in
            ASSERT( n_positive == 0 ||  n_positive == 3);
            return { s3_dim_starts[0]+i_vertex, s3_dim_starts[0]+i_vertex};
        }
    }
}


std::vector<std::vector<arma::uvec>> ComputeIntersection<2,3>::_on_faces()
{
    std::vector<std::vector<arma::uvec>> on_faces;
    
    on_faces.resize(3);
    arma::uvec::fixed<RefElement<3>::n_sides> v; v.zeros();
    
    on_faces[0].resize(RefElement<3>::n_nodes, v);
    for(uint i=0; i<RefElement<3>::n_nodes; i++)
        for(uint j=0; j<RefElement<3>::n_sides_per_node; j++)
            on_faces[0][i](RefElement<3>::interact(Interaction<2,0>(i))[j]) = 1;
        
    on_faces[1].resize(RefElement<3>::n_lines, v);
    for(uint i=0; i<RefElement<3>::n_lines; i++)
        for(uint j=0; j<RefElement<3>::n_sides_per_line; j++)
            on_faces[1][i](RefElement<3>::interact(Interaction<2,1>(i))[j]) = 1;
    
    on_faces[2].resize(RefElement<3>::n_sides, v);
    for(uint i=0; i<RefElement<3>::n_sides; i++)
        on_faces[2][i](i) = 1;
            
//     DBGCOUT("Print on_faces:\n");
//     for(uint d=0; d<on_faces.size(); d++)
//         for(uint i=0; i<on_faces[d].size(); i++){
//             for(uint j=0; j<on_faces[d][i].n_elem; j++)
//                 cout << on_faces[d][i](j) << " ";
//             cout << endl;
//         }
    return on_faces;
}


void ComputeIntersection<2,3>::print_plucker_coordinates(std::ostream &os){
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
}

void ComputeIntersection<2,3>::print_plucker_coordinates_tree(std::ostream &os){
	os << "ComputeIntersection<2,3> Plucker Coordinates Tree:" << endl;
	print_plucker_coordinates(os);
	for(unsigned int i = 0; i < 6;i++){
		os << "ComputeIntersection<1,2>["<< i <<"] Plucker Coordinates:" << endl;
		CI12[i].print_plucker_coordinates(os);
	}
	for(unsigned int i = 0; i < 3;i++){
		CI13[i].print_plucker_coordinates_tree(os);
	}
}

