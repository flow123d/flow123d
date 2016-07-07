/*!
 *
﻿* Copyright (C) 2015 Technical University of Liberec.  All rights reserved.
 * 
 * This program is free software; you can redistribute it and/or modify it under
 * the terms of the GNU General Public License version 3 as published by the
 * Free Software Foundation. (http://www.gnu.org/licenses/gpl-3.0.en.html)
 * 
 * This program is distributed in the hope that it will be useful, but WITHOUT
 * ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 * FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.
 *
 * 
 * @file    plucker.hh
 * @brief   Plucker coordinates class.
 * @author  Viktor Fris, Pavel Exner
 *
 */

#include <armadillo>
#include <iostream>
#include "system/system.hh"

namespace computeintersection{

#ifndef _PLUCKER_H
#define _PLUCKER_H

/** @brief Plucker coordinates representing line given by points A,B.
 * 
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
 * question: is compute() used from outside ? 
 * answer: Yes
 */
class Plucker{
private:

	arma::vec6 coordinates_; ///< Plucker coordinates.
	bool computed_;          ///< True, if Plucker coordinates are computed; false otherwise.

public:
    /** Default constructor.
     * Plucker coordinates are not computed.
     */
	Plucker();
	/** @brief Creates Plucker coordinates of a line AB.
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

    /** @brief Compute Plucker coordinates and set computed to true.
     * @param a A point from AB line
     * @param b B point from AB line (order of points determines orientation!)
     */
	void compute(const arma::vec3 &a, const arma::vec3 &b);

	/// Gets directional vector U.
	arma::vec3 get_u_vector() const;

	/// Gets cross product vector UxA.
	arma::vec3 get_ua_vector() const;

    /// Gets Plucker coordinates.
	arma::vec6 get_plucker_coords() const;
    
    /// Friend output operator.
    friend std::ostream& operator<< (std::ostream& os, const Plucker& p);
};

/// Operator for printing Plucker coordinates.
std::ostream& operator<<(std::ostream& os, const Plucker& p);

/****************** inline implementation *****************************/
inline Plucker::Plucker()
{   computed_ = false; }

inline Plucker::Plucker(const arma::vec3 &a,const arma::vec3 &b)
{   compute(a, b);
    computed_ = true; }

inline Plucker::Plucker(const Plucker &p)
{   coordinates_ = p.get_plucker_coords();
    computed_ = p.is_computed(); }

inline double Plucker::operator[](const unsigned int index) const
{   ASSERT_DBG(computed_);
    return coordinates_[index]; }

inline void Plucker::clear()
{   computed_ = false; }

inline bool Plucker::is_computed() const
{   return computed_; }

inline arma::vec3 Plucker::get_u_vector() const
{   ASSERT_DBG(computed_);
    return coordinates_(arma::span(0,2)); }

inline arma::vec3 Plucker::get_ua_vector() const
{   ASSERT_DBG(computed_);
    return coordinates_(arma::span(3,5)); }

inline arma::vec6 Plucker::get_plucker_coords() const
{   ASSERT_DBG(computed_);
    return coordinates_; }

#endif

} // END namespace_close
