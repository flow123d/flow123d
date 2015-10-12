#include <armadillo>
#include <iostream>

using namespace std;
namespace computeintersection{

#ifndef _PLUCKER_H
#define _PLUCKER_H

/**
 * Plucker class represents a line by 6 dimensional vector.
 * After inserting a three-dimensional points A and B, which represents the line,
 * class creates plucker coordinates of the line.
 *
 * Class also can compute a product of two plucker coordinates.
 * 
 * Description of Plücker Coordinates:
 * https://en.wikipedia.org/wiki/Pl%C3%BCcker_coordinates
 *
 * Empty constructor is used for passing object to pointers from different places
 * coordinates data are filled after calling method "compute"
 * a flag "computed" is for comparison if coordinates data are filled
 *
 * TODO: 
 * - is compute() used from outside ? Yes
 */
class Plucker{
private:

	arma::vec6 coordinates;
	bool computed;

public:
	Plucker();
	Plucker(const arma::vec3 &a, const arma::vec3 &b);
	/// copy constructor
	Plucker(const Plucker &p);
	inline ~Plucker(){};

	inline double operator[](const unsigned int index) const{
		return coordinates[index];
	};

	/// Compute product of two Plücker coordinates
	double operator*(const Plucker &b);

	inline void clear(){computed = false;};

	inline bool is_computed() const{
		return computed;
	};


	// Compute Plücker coordinates and set computed to true
	void compute(const arma::vec3 &a, const arma::vec3 &b);

	/// get directional vector U
	inline arma::vec3 get_u_vector() const{
		return coordinates(arma::span(0,2));
	};

	/// get cross product vector UxA
	inline arma::vec3 get_ua_vector() const{
		return coordinates(arma::span(3,5));
	};

	inline arma::vec6 get_plucker_coords() const{
		return coordinates;
	};

    ///TODO: convert to operator <<
	void toString();
};

#endif

} // END namespace_close
