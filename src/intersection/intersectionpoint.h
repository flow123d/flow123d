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

	int edge_idx;
	int side_idx;

	public:
	inline IntersectionPoint(const arma::vec::fixed<N+1> &lc1,
					  const arma::vec::fixed<M+1> &lc2,
					  int edge = -1,
					  int side = -1)
					  : local_coords1(lc1),
					    local_coords2(lc2),
					    edge_idx(edge),
					    side_idx(side){};
	inline arma::vec::fixed<N+1> getLocalCoords1(){
			return local_coords1;
		};
	inline arma::vec::fixed<M+1> getLocalCoords2(){
			return local_coords2;
		};

	inline void setEdge(int e){
		edge_idx = e;
	};

	inline void setSide(int s){
		side_idx = s;
	};

	inline int getEdge(){
		return edge_idx;
	};

	inline int getSide(){
		return side_idx;
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
