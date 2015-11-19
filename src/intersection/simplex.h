#ifndef SIMPLEX_H_
#define SIMPLEX_H_

#include <armadillo>
#include "system/system.hh"

namespace computeintersection {

/**
 * Simplex<N> represents N-dimensional simplex,
 * Simplex<0> = pointer to 3D point
 * Simplex<1> = abscissa
 * Simplex<2> = triangle
 * Simplex<3> = tetrahedron
 * Sub - simplices are made in lexicographical order
 *
 * 				Simplex<3> with 4 points 0,1,2,3 creates:
 * 								|
 * 					 -----------------------------------------------------------------------------------------------------------------------------------------------
 * 					|												|												|												|
 * 					S<2>(0,1,2) 									S<2>(0,1,3) 									S<2>(0,2,3) 									S<2>(1,2,3)
 * 						|												|												|												|
 * 		 -------------------------------				 -------------------------------				 -------------------------------				 -------------------------------
 * 		|				|				|				|				|				|				|				|				|				|				|				|
 * 		S<1>(0,1) 		S<1>(0,2) 		S<1>(1,2) 		S<1>(0,1) 		S<1>(0,3) 		S<1>(1,3) 		S<1>(0,2) 		S<1>(0,3) 		S<1>(2,3) 		S<1>(1,2) 		S<1>(1,3) 		S<1>(2,3)
 * 			|				|				|				|				|				|				|				|				|				|				|				|
 * 		 -------		 -------		 -------		 -------		 -------		 -------		 -------		 -------		 -------		 -------		 -------		 -------
 * 		|		|		|		|		|		|		|		|		|		|		|		|		|		|		|		|		|		|		|		|		|		|		|		|
 *		S<0>(0) S<0>(1) S<0>(0) S<0>(2) S<0>(1) S<0>(2) S<0>(0) S<0>(1) S<0>(0) S<0>(3) S<0>(1) S<0>(3) S<0>(0) S<0>(2) S<0>(0) S<0>(3) S<0>(2) S<0>(3) S<0>(1) S<0>(2) S<0>(1) S<0>(3) S<0>(2) S<0>(3)
 *
 * Simplex<0> has pointer to 3D point, because 3D point will be only once in memory
 *
 * https://en.wikipedia.org/wiki/Simplex
 */
template<unsigned int N> class Simplex;

/// Operator for printing Simplex<N>.
template<unsigned int N> std::ostream& operator<<(std::ostream& os, const Simplex<N>& s);

/**
 * Simplex<0> represents a point in 3D
 * it has a pointer to coordinates of point in 3D because of mesh implementation
 */
template<> class Simplex<0> {
private:
	arma::vec3* coords_;    ///< Point coordinates.
public:
	Simplex();  /// Default constructor. Does not set coordinates pointer.
	
    /** Constuctor that sets the point coordinates.
	 * @param field - array of pointers to point coordinates;
	 * it takes just the first element of the input array for case of Simplex<0>
	 */
	Simplex(arma::vec3 **field);
    
	~Simplex(); ///< Destructor.

	/**
	 * Setter for point coordinates
	 */
	void set_simplices(arma::vec3 **field);

    /// Returns the point coordinates.
	arma::vec3 &point_coordinates();

	/// Friend output operator.
	friend std::ostream& operator<< <>(std::ostream& os, const Simplex<0>& s);
};

template<unsigned int N> class Simplex {
private:
    /**
     * Every simplex<N> has (N+1) simplices of dimension (N-1)
     * TODO: (idea - when used in new mesh) replace with references 
     */
	Simplex<N - 1> simplices_[N + 1];
public:
	Simplex();  ///< Default (empty) constructor.

	/** Constuctor that sets the coordinates of all vertices.
	 * @param field - array of pointers to point coordinates of the vertices
	 * array must contain (N+1) elements for case Simplex<N>
	 */
	Simplex(arma::vec3 **field_of_pointers_to_coordinates);

    ~Simplex();    ///< Destructor.
    
	/// Creating sub-simplices in lexicografic order
	void set_simplices(arma::vec3 **field_of_pointers_to_coordinates);

    /// Returns subsimplex of index @p idx.
	Simplex<N - 1> &operator[](unsigned int idx);

	/// Get simplex of abscissa from different simplices - if it has own implementation in .cpp file
	Simplex<1> &abscissa(unsigned int idx);
    
    /// Friend output operator.
    friend std::ostream& operator<< <>(std::ostream& os, const Simplex<N>& s);
};



/********************************************* IMPLEMENTATION ***********************************************/

inline Simplex< 0 >::Simplex()
{   coords_ = nullptr; }

inline Simplex< 0  >::Simplex(arma::vec3** field)
{   ASSERT(field != nullptr, "Null pointer given in the constructor.");
    coords_ = field[0]; }

inline Simplex< 0  >::~Simplex(){}

inline void Simplex< 0  >::set_simplices(arma::vec3** field)
{   ASSERT(field != nullptr, "Null pointer given in the setter.");
    coords_ = field[0]; }

inline arma::vec3& Simplex< 0  >::point_coordinates()
{   ASSERT(coords_ != nullptr, "Null pointer given in the constructor.");
    return *coords_; }



template<unsigned int N> Simplex<N>::Simplex(){}

template<unsigned int N> Simplex<N>::Simplex(arma::vec3 **field_of_pointers_to_coordinates)
{   set_simplices(field_of_pointers_to_coordinates); }

template<unsigned int N> Simplex<N>::~Simplex(){}

template<unsigned int N> Simplex<N-1>& Simplex<N>::operator[](unsigned int idx)
{   ASSERT(idx < N+1, "Index out of bounds (number of subsimplices.)");
    return simplices_[idx]; }

} // END namespace_close
#endif /* SIMPLEX_H_ */
