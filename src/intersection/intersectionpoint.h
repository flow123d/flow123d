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

	public:

	inline IntersectionPoint(){
		side_idx1 = -1;
		side_idx2 = -1;
	};
	inline IntersectionPoint(const arma::vec::fixed<N+1> &lc1,
					  const arma::vec::fixed<M+1> &lc2,
					  int side1 = -1,
					  int side2 = -1)
					  : local_coords1(lc1),
					    local_coords2(lc2),
					    side_idx1(side1),
					    side_idx2(side2){};
	inline ~IntersectionPoint(){};

	inline void clear(){
		local_coords1.zeros();
		local_coords2.zeros();
		side_idx1 = -1;
		side_idx2 = -1;
	};

	inline void print(){
		cout << "Local coords on the first element on side(" << side_idx1 << ")" << endl;
		local_coords1.print();
		cout << "Local coords on the second element on side(" << side_idx2 << ")" << endl;
		local_coords2.print();
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

	inline int getSide1(){
		return side_idx1;
	};

	inline int getSide2(){
		return side_idx2;
	};

	inline int getDim1(){
		return N;
	};
	inline int getDim2(){
		return M;
	};

};


/*
template<> class IntersectionPoint<2,1>{

	arma::vec3 local_coords1;
	double local_coords2;
	unsigned int index_edge;

	public:
		IntersectionPoint(const arma::vec3 &lc1, const double &lc2, unsigned int ie) : local_coords1(lc1), local_coords2(lc2),index_edge(ie){};
	inline arma::vec3 &local_coords_triangle(){return local_coords1;};
	inline double &local_coords_triangle(unsigned int index){return local_coords1[index];};
	inline double &local_coord_abscissa(){return local_coords2;};
	inline unsigned int &edge(){return index_edge;};
};

template<> class IntersectionPoint<1,2>{

	double local_coords1;
	arma::vec3 local_coords2;

	unsigned int index_edge;
	unsigned int index_side;

public:
	IntersectionPoint(const double &lc1,const arma::vec3 &lc2, unsigned int ie, unsigned int is) : local_coords1(lc1), local_coords2(lc2), index_edge(ie), index_side(is){};
	inline double &local_coords_abscissa(){return local_coords1;};
	inline arma::vec3 &local_coords_triangle(){return local_coords2;};
	inline double &local_coord_triangle(unsigned int index){return local_coords2[index];};
	inline unsigned int &edge(){return index_edge;};
	inline unsigned int &side(){return index_side;};
};

template<> class IntersectionPoint<1,3>{

	double local_coords1;
	arma::vec3 local_coords2;

	unsigned int index_edge;

public:
	IntersectionPoint(const double &lc1,const arma::vec3 &lc2, unsigned int ie) : local_coords1(lc1), local_coords2(lc2), index_edge(ie){};
	inline double &local_coords_abscissa(){return local_coords1;};
	inline arma::vec3 &local_coords_tetraheadron(){return local_coords2;};
	inline double &local_coord_tetraheadron(unsigned int index){return local_coords2[index];};
	inline unsigned int &edge(){return index_edge;};
};

template<> class IntersectionPoint<2,3>{

	arma::vec3 local_coords1;
	arma::vec3 local_coords2;

public:
	IntersectionPoint(const arma::vec3 &lc1,const arma::vec3 &lc2) : local_coords1(lc1), local_coords2(lc2){};
};*/


//template<int M, int N> arma::vec::fixed<N+1> interpolate(arma::vec::fixed<M+1> &coord, unsigned int sub_simplex_idx);
// interpolation of barycentric coordinate on Simplex<M> to Simplex<N>

} // END namespace
