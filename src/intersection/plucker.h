#include <armadillo>
#include <iostream>

using namespace std;
namespace computeintersection{

template<int N> class Simplex;

class Plucker{
private:
	arma::vec6 coordinates;
	//double coordinates[6];
public:
	inline Plucker(){};
	//Plucker(arma::vec3 u, arma::vec3 a);
	Plucker(const arma::vec3 u, const arma::vec3 a);
	inline ~Plucker(){};
	inline double operator[](int index){return coordinates[index];};
	double operator*(Plucker b);
	void operator*(double number);
	void toString();
};
} // END namespace_close
