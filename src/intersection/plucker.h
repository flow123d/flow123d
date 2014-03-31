#include <armadillo>
#include <iostream>

using namespace std;
namespace computeintersection{

class Plucker{
private:
	arma::vec6 coordinates;
	bool computed;

public:
	Plucker();
	Plucker(const arma::vec3 u, const arma::vec3 a);
	inline ~Plucker(){};

	inline double operator[](int index){return coordinates[index];};

	double operator*(Plucker b);
	void operator*(double number);

	void setComputed(bool computed);
	void clear();
	bool isComputed();


	void compute(const arma::vec3 u, const arma::vec3 a);

	void toString();
};
} // END namespace_close
