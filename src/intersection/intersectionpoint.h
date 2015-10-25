#ifndef INTERSECTIONPOINT_H_
#define INTERSECTIONPOINT_H_

#include <armadillo>
#include <array>
#include "system/system.hh"
#include "refsimplex.h"

using namespace std;
namespace computeintersection{

/**
 * Have separate class template for IntersectionPoint as it appears at output
 * and other internal class template e.g. TriangleLineIntersection for internal
 * intersections with additional info.
 *
 * 
 * Represents a point from simplex<N> and simplex<M> intersection
 * contains bary coords of a point on simplex<N> and simplex<M>
 *
 */
template<int N, int M> class IntersectionPoint {

	arma::vec::fixed<N+1> local_coords1; // bary coords of a point on simplex<N>
	arma::vec::fixed<M+1> local_coords2; // bary coords of a point on simplex<M>

	int side_idx1; // For case N = 2, M = 3 -> index of a triangle line
	int side_idx2; // For case N = 1 or N = 2, M = 3 -> index of a tetrahedron side

	unsigned int orientation; // orientation from intersection using plucker coords

	bool is_vertex_; // point is a vertex of triangle
	bool is_patological_; // points is a vertex of tetrahedron or it is in side of tetrahedron or it is in edge of tetrahedron

	public:

	inline IntersectionPoint(){
		clear();
	};
	inline IntersectionPoint(const arma::vec::fixed<N+1> &lc1,
					  const arma::vec::fixed<M+1> &lc2,
					  int side1 = -1,
					  int side2 = -1,
					  unsigned int ori = 1,
					  bool vertex = false,
					  bool patological = false)
					  : local_coords1(lc1),
					    local_coords2(lc2),
					    side_idx1(side1),
					    side_idx2(side2),
					    orientation(ori),
					    is_vertex_(vertex),
					    is_patological_(patological){};
	inline ~IntersectionPoint(){};

	/// Fliping dimension of an intersectionpoint
	inline IntersectionPoint(IntersectionPoint<M, N> &IP){
		local_coords1 = IP.get_local_coords2();
		local_coords2 = IP.get_local_coords1();
		side_idx1 = IP.get_side2();
		side_idx2 = IP.get_side1();
		orientation = IP.get_orientation();
		is_vertex_ = IP.is_vertex();
		is_patological_ = IP.is_patological();
	};

	/* Constructor interpolates the second bary coords of IntersectionPoint<N,M-1> to IntersectionPoint<N,M>
     * Allowed only from dimension 1 to 2 and from 2 to 3
     * */
	inline IntersectionPoint(IntersectionPoint<N,M-1> &IP){
		arma::vec::fixed<M+1> interpolated;
		//TODO: PE; try to replace if cases; interpolate<M-1>
		if(M == 3){
			interpolated = RefSimplex<3>::interpolate<2>(IP.get_local_coords2(), IP.get_side2());
		}else if(M == 2){
			interpolated = RefSimplex<2>::interpolate<1>(IP.get_local_coords2(), IP.get_side2());
		}else{
			ASSERT(false,"Wrong the second dimension in an IntersectionPoint (allowed 2 and 3 only)");
			interpolated.zeros();
		}
		local_coords1 = IP.get_local_coords1();
		local_coords2 = interpolated;
		side_idx1 = IP.get_side1();
		side_idx2 = IP.get_side2();
		orientation = IP.get_orientation();
		is_vertex_ = IP.is_vertex();
		is_patological_ = IP.is_patological();
	}

	/* Constructor interpolates the second bary coords of IntersectionPoint<N,M-2> to IntersectionPoint<N,M>
	 * Allowed only from dimension 1 to 3
	 * */
	inline IntersectionPoint(IntersectionPoint<N,M-2> &IP){
     	arma::vec::fixed<M+1> interpolated;
     	if(M == 3){
     		interpolated = RefSimplex<3>::interpolate<1>(IP.get_local_coords2(), IP.get_side2());
     	}else{
     		ASSERT(false,"Wrong the second dimension in an IntersectionPoint (allowed 3 only)");
     		interpolated.zeros();
     	}
     	local_coords1 = IP.get_local_coords1();
		local_coords2 = interpolated;
		side_idx1 = IP.get_side1();
		side_idx2 = IP.get_side2();
		orientation = IP.get_orientation();
		is_vertex_ = IP.is_vertex();
		is_patological_ = IP.is_patological();
	}

	inline void clear(){
		local_coords1.zeros();
		local_coords2.zeros();
		side_idx1 = -1;
		side_idx2 = -1;
		orientation = 1;
		is_vertex_ = false;
		is_patological_ = false;
	};

	inline void print(){
		cout << "Local coords on the first element on side(" << side_idx1 << ")" << endl;
		local_coords1.print();
		cout << "Local coords on the second element on side(" << side_idx2 << ")" << endl;
		local_coords2.print();
		cout << "Orientation: " << orientation << " Vertex: " << is_vertex_ << " Patological: " << is_patological_ << endl;
	};

	inline const arma::vec::fixed<N+1> &get_local_coords1() const{
		return local_coords1;
	};
	inline const arma::vec::fixed<M+1> &get_local_coords2() const{
		return local_coords2;
	};

	inline void set_side1(int s){
		side_idx1 = s;
	};

	inline void set_side2(int s){
		side_idx2 = s;
	};

	inline void set_orientation(unsigned int o){
		orientation = o;
	};

	inline void set_is_vertex(bool iv){
		is_vertex_ = iv;
	}

	inline void set_is_patological(bool ip){
		is_patological_ = ip;
	}

	inline int get_side1() const{
		return side_idx1;
	};

	inline int get_side2() const{
		return side_idx2;
	};

	inline unsigned int get_orientation() const{
		return orientation;
	};

	inline bool is_vertex() const{
		return is_vertex_;
	};

	inline bool is_patological() const{
		return is_patological_;
	};

	inline friend ostream& operator<<(ostream& os, const IntersectionPoint<N,M>& IP){
		os << "test";
		return os;
	};
	/**
	 * For convex hull polygon tracing
	 */
	bool operator<(const IntersectionPoint<N,M> &ip) const;
};


} // END namespace
#endif /* INTERSECTIONPOINT_H_ */
