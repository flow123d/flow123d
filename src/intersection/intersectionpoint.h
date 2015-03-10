#include <armadillo>
#include <array>
#include "system/system.hh"
#include "refsimplex.h"

using namespace std;
namespace computeintersection{

//class ProlongationPoint;

/**
 * Class doc.
 * Naming convention.
 *
 * Have separate class template for IntersectionPoint as it appears at output
 * and other internal class template e.g. TriangleLineIntersection for internal
 * intersections with additional info.
 *
 * Describe data attributes.
 */
template<int N, int M> class IntersectionPoint {

	arma::vec::fixed<N+1> local_coords1;
	arma::vec::fixed<M+1> local_coords2;

	int side_idx1;
	int side_idx2;

	unsigned int orientation;

	bool is_vertex;
	bool is_patological;

	public:


	// Possibly just call clear()
	inline IntersectionPoint(){
		this->clear();
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
					    is_vertex(vertex),
					    is_patological(patological){};
	inline ~IntersectionPoint(){};

	inline void clear(){
		local_coords1.zeros();
		local_coords2.zeros();
		side_idx1 = -1;
		side_idx2 = -1;
		orientation = 1;
		is_vertex = false;
		is_patological = false;
	};

	inline void print(){
		cout << "Local coords on the first element on side(" << side_idx1 << ")" << endl;
		local_coords1.print();
		cout << "Local coords on the second element on side(" << side_idx2 << ")" << endl;
		local_coords2.print();
		cout << "Orientation: " << orientation << " Vertex: " << is_vertex << " Patological: " << is_patological << endl;
	};

	inline double get_lc1_coord(unsigned int i) const{
		return local_coords1[i];
	};

	inline arma::vec::fixed<N+1> &get_local_coords1(){
		return local_coords1;
	};
	inline arma::vec::fixed<M+1> &get_local_coords2(){
		return local_coords2;
	};

	inline void setSide1(int s){
		side_idx1 = s;
	};

	inline void setSide2(int s){
		side_idx2 = s;
	};

	inline void setOrientation(unsigned int o){
		orientation = o;
	};

	inline void setIsVertex(bool iv){
		is_vertex = iv;
	}

	inline void setIsPatological(bool ip){
		is_patological = ip;
	}

	inline int getSide1() const{
		return side_idx1;
	};

	inline int getSide2() const{
		return side_idx2;
	};

	inline unsigned int getOrientation() const{
		return orientation;
	};

	inline bool isVertex() const{
		return is_vertex;
	};

	inline bool isPatological() const{
		return is_patological;
	};

	bool operator<(const IntersectionPoint<N,M> &ip) const;
};

class TriangleLineIntersections{

	int side_idx1;
	int side_idx2;

	unsigned int orientation;

	bool is_vertex;
	bool is_patological;
};

} // END namespace
