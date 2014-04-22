#include <armadillo>
#include "system/system.hh"

using namespace std;
namespace computeintersection{

//class ProlongationPoint;

template<int N, int M> class IntersectionPoint{};

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
};

} // END namespace
