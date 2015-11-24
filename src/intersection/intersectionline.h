/*
 * IntersectionLine.h
 *
 *  Created on: 8.4.2015
 *      Author: viktor
 */

#ifndef INTERSECTIONLINE_H_
#define INTERSECTIONLINE_H_

#include "system/system.hh"

namespace computeintersection{

//forwward declare
template<unsigned int, unsigned int> class IntersectionPoint;

/**
 * IntersectionLine represents an intersection of a line(Simplex<1>) and a tetrahedron(Simplex<3>).
 * Contains array of intersection points (necessary not in specified order before calling @p trace_line).
 * TODO: what about keeping directly ElementFullIterator ?
 */
class IntersectionLine{

	std::vector<IntersectionPoint<1,3>> i_points_;
	unsigned int element_1d_idx;
	unsigned int element_3d_idx;

public:

	IntersectionLine(){};                                       ///< Default constructor.
	IntersectionLine(unsigned int ele_1D, unsigned int ele_3D)  ///< Constructor taking in element indices.
        :element_1d_idx(ele_1D), element_3d_idx(ele_3D){};
	~IntersectionLine(){};                                      ///< Destructor.

	/// Returns intersection points by reference.
	std::vector<IntersectionPoint<1,3>> &points();

	/// Returns intersection points by constant reference.
	const std::vector<IntersectionPoint<1,3>> &points() const;

	/**
	 * Switches intersection points if they are not in correct order.
	 */
// 	void trace_line();

    /// Returns intersection point of given @p index.
    const IntersectionPoint<1,3> &operator[](unsigned int index) const;

    
    unsigned int size() const;          ///< Returns number of intersection points.
    unsigned int ele_1d_idx() const;    ///< Returns index of 1D element.
	unsigned int ele_3d_idx() const;    ///< Returns index of 3D element.
    
    /// Computes the relative length of intersection line.
    double compute_length();

};

/********************************************* IMPLEMENTATION ***********************************************/

inline std::vector< IntersectionPoint< 1, 3 > >& IntersectionLine::points()
{   return i_points_; }

inline const std::vector< IntersectionPoint< 1, 3 > >& IntersectionLine::points() const
{   return i_points_; }

inline const IntersectionPoint< 1, 3 >& IntersectionLine::operator[](unsigned int index) const
{   ASSERT(index < i_points_.size(), "Index out of bounds.");
    return i_points_[index]; }

inline unsigned int IntersectionLine::size() const
{   return i_points_.size(); }

inline unsigned int IntersectionLine::ele_1d_idx() const
{   return element_1d_idx; }

inline unsigned int IntersectionLine::ele_3d_idx() const
{   return element_3d_idx; }


} // END NAMESPACE
#endif /* INTERSECTIONLINE_H_ */
