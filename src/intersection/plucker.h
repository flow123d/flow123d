#include <armadillo>
// #include <iostream>

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

	arma::vec6 coordinates; ///< Plucker coordinates.
	bool computed;          ///< True, if Plucker coordinates are computed; false otherwise.

public:
    /** Default constructor.
     * Plucker coordinates are not computed.
     */
	Plucker();
	/**
     * Creates Plucker coordinates of a line AB.
	 * @param a - A point from AB line
	 * @param b - B point from AB line
	 */
	Plucker(const arma::vec3 &a, const arma::vec3 &b);
	/// Copy constructor.
	Plucker(const Plucker &p);
    /// Destructor.
	~Plucker(){};

    /// Returns Plucker coordinate of @p index.
	double operator[](const unsigned int index) const;

	/// Compute product of two Plucker coordinates.
	double operator*(const Plucker &b);

    /// Sets the flag computed on false.
	void clear();

    /// Return true if Plucker coordinates have been computed already.
	bool is_computed() const;

	/// Compute Plucker coordinates and set computed to true.
	void compute(const arma::vec3 &a, const arma::vec3 &b);

	/// Gets directional vector U.
	arma::vec3 get_u_vector() const;

	/// Gets cross product vector UxA.
	arma::vec3 get_ua_vector() const;

    /// Gets Plucker coordinates.
	arma::vec6 get_plucker_coords() const;

    ///TODO: convert to operator <<
	void toString();
};



/****************** inline implementation *****************************/
inline double Plucker::operator[](const unsigned int index) const
{ return coordinates[index]; }


inline void Plucker::clear()
{ computed = false; }

inline bool Plucker::is_computed() const
{ return computed; }

inline arma::vec3 Plucker::get_u_vector() const
{ return coordinates(arma::span(0,2)); }

inline arma::vec3 Plucker::get_ua_vector() const
{ return coordinates(arma::span(3,5)); }

inline arma::vec6 Plucker::get_plucker_coords() const
{ return coordinates; }

#endif

} // END namespace_close
