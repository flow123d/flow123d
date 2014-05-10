#include <armadillo>
#include <array>
#include "system/system.hh"
#include "refsimplex.h"

using namespace std;
namespace computeintersection{

//class ProlongationPoint;


template<int N, int M> class IntersectionPoint{

	arma::vec::fixed<N+1> local_coords1;
	arma::vec::fixed<M+1> local_coords2;

	int side_idx1;
	int side_idx2;

	unsigned int orientation;
	bool is_vertex;

	public:

	inline IntersectionPoint(){
		side_idx1 = -1;
		side_idx2 = -1;
		orientation = 1;
		is_vertex = false;
	};
	inline IntersectionPoint(const arma::vec::fixed<N+1> &lc1,
					  const arma::vec::fixed<M+1> &lc2,
					  int side1 = -1,
					  int side2 = -1,
					  unsigned int ori = 1,
					  bool vertex = false)
					  : local_coords1(lc1),
					    local_coords2(lc2),
					    side_idx1(side1),
					    side_idx2(side2),
					    orientation(ori),
					    is_vertex(vertex){};
	inline ~IntersectionPoint(){};

	inline void clear(){
		local_coords1.zeros();
		local_coords2.zeros();
		side_idx1 = -1;
		side_idx2 = -1;
		orientation = 1;
		is_vertex = false;
	};

	inline void print(){
		cout << "Local coords on the first element on side(" << side_idx1 << ")" << endl;
		local_coords1.print();
		cout << "Local coords on the second element on side(" << side_idx2 << ")" << endl;
		local_coords2.print();
		cout << "Orientation: " << orientation << endl;
	};

	inline arma::vec::fixed<N+1> getLocalCoords1(){
			return local_coords1;
		};
	inline arma::vec::fixed<M+1> getLocalCoords2(){
			return local_coords2;
		};

	inline void setLocalCoords1(const arma::vec::fixed<N+1> lc1){
			local_coords1 = lc1;
	};

	inline void setLocalCoords2(const arma::vec::fixed<M+1> lc2){
			local_coords2 = lc2;
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

	inline int getSide1(){
		return side_idx1;
	};

	inline int getSide2(){
		return side_idx2;
	};

	inline unsigned int getOrientation(){
		return orientation;
	};

	inline bool isVertex(){
		return is_vertex;
	};

	inline int getDim1(){
		return N;
	};
	inline int getDim2(){
		return M;
	};

};

} // END namespace
