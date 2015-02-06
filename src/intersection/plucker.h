#include <armadillo>
#include <iostream>

using namespace std;
namespace computeintersection{

#ifndef _PLUCKER_H
#define _PLUCKER_H

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
	bool isComputed() const;


	void compute(const arma::vec3 u, const arma::vec3 a);

	// get directional vector U
	arma::vec3 getU();

	// get cross product vector UxA
	arma::vec3 getUA();

	arma::vec6 getPlucker() const;

	void setPlucker(const Plucker &p);

	void toString();
};

#endif

} // END namespace_close
