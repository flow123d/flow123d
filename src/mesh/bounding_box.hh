/*!
 *
 * Copyright (C) 2007 Technical University of Liberec.  All rights reserved.
 *
 * Please make a following refer to Flow123d on your project site if you use the program for any purpose,
 * especially for academic research:
 * Flow123d, Research Centre: Advanced Remedial Technologies, Technical University of Liberec, Czech Republic
 *
 * This program is free software; you can redistribute it and/or modify it under the terms
 * of the GNU General Public License version 3 as published by the Free Software Foundation.
 *
 * This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
 * without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 * See the GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License along with this program; if not,
 * write to the Free Software Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA 021110-1307, USA.
 *
 *
 * $Id: bounding_box.hh 1567 2012-02-28 13:24:58Z jan.brezina $
 * $Revision: 1567 $
 * $LastChangedBy: jan.brezina $
 * $LastChangedDate: 2012-02-28 14:24:58 +0100 (Tue, 28 Feb 2012) $
 *
 *
 */

#ifndef BOX_ELEMENT_HH_
#define BOX_ELEMENT_HH_

#include "system/system.hh"
#include "system/global_defs.h"
#include "mesh/point.hh"

#include <vector>
#include <armadillo>

/**
 * @brief Bounding box in 3d ambient space.
 *
 * Primary intention is usage in BIHTree and various speedups of non-compatible
 * intersections.
 *
 * Copy constructor and assignment are default provided by compiler.
 * These can be used to set bounds latter on without particular method
 * to this end:
 *
 * @code
 * BoundingBox box;			// non-initialized box
 * box=BoundingBox( arma::vec3("0 1 2"), arma::vec3("4 5 6") );
 * @endcode
 *
 * Don;t worry about performance, all is inlined.
 */
class BoundingBox {
public:
	TYPEDEF_ERR_INFO( EI_split_point, double);
	TYPEDEF_ERR_INFO( EI_interval_left, double);
	TYPEDEF_ERR_INFO( EI_interval_right, double);
	DECLARE_EXCEPTION(ExcSplitting, << "Split point " << EI_split_point::val
			                        << "out of bounds: <" << EI_interval_left::val
			                        << ", " << EI_interval_right::val << ">\n");

	/// Currently we set dimension to 3.
	static const unsigned int dimension = 3;
	/// Currently we assume
	typedef Space<dimension>::Point Point;

	/**
	 * Default constructor.
	 * No initialization of vertices. Be very careful using this.
	 * One necessary usage is vector of BoundigBox.
	 */
	BoundingBox() {}

	/**
	 * Constructor for point box.
	 */
	BoundingBox(const Point &min)
	: min_vertex_(min), max_vertex_(min)
	{};

	/**
	 * Constructor.
	 *
	 * From given minimal and maximal vertex.
	 */
	BoundingBox(const Point &min, const Point &max)
	: min_vertex_(min), max_vertex_(max)
	{
		ASSERT( arma::min( min <= max ) , "Wrong coordinates in constructor.");
	};

	/**
	 * Constructor.
	 *
	 * Make bounding box for set of points.
	 */
	BoundingBox(const vector<Point> &points);

	/**
	 * Set maximum in given axis.
	 */
	void set_max(unsigned int axis, double max) {
		ASSERT_LESS(axis , dimension);
		ASSERT_LE( min(axis) , max);
		max_vertex_[axis] = max;
	}

	/**
	 * Set minimum on given axis.
	 */
	void set_min(unsigned int axis, double min) {
		ASSERT_LESS(axis, dimension);
		ASSERT_LE(min , max(axis));
		min_vertex_[axis] = min;
	}

	/**
	 *  Return minimal vertex of the bounding box.
	 */
	const Point &min() const {
		return min_vertex_;
	}

	/**
	 *  Return maximal vertex of the bounding box.
	 */
	const Point &max() const {
		return max_vertex_;
	}

	/**
	 *  Return minimal value on given axis.
	 */
	double min(unsigned int axis) const {
		return min()[axis];
	}

	/**
	 *  Return maximal value on given axis.
	 */
	double max(unsigned int axis) const {
		return max()[axis];
	}


	/**
	 *  Return size of the box in given axis.
	 */
	double size(unsigned int axis) const {
		return max()[axis] - min()[axis];
	}


	/**
	 *  Return center of the bounding box.
	 */
	Point center() const {
		return (max_vertex_ + min_vertex_) / 2.0;
	}

     /**
     * Return center of projection of the bounding box to given @p axis.
     * Axis coding is: 0 - axis x, 1 - axis y, 2 - axis z.
     */
    double projection_center(unsigned int axis) const {
    	ASSERT_LESS(axis, dimension);
    	return (max_vertex_[axis] + min_vertex_[axis])/2;
    }

    /**
     * Returns true is the  box element contains @p point
     *
     * @param point Testing point
     * @return True if box element contains point
     */
    bool contains_point(const Point &point) const
    {
    	for (unsigned int i=0; i<dimension; i++) {
    		if ((point(i) + epsilon < min_vertex_(i)) ||
    			(point(i) > epsilon + max_vertex_(i))) return false;
    	}
    	return true;
    }

    /**
     * Returns true if two bounding boxes have intersection.
     *
     * This serves as an estimate of intersection of elements.
     * To make it safe (do not exclude possible intersection) for
     * 1d and 2d elements aligned with axes, we use some tolerance.
     * Since this tolerance is fixed, there could be problem with
     * highly refined meshes (get false positive result).
     */
    bool intersect(const BoundingBox &b2) const
    {
    	for (unsigned int i=0; i<dimension; i++) {
    		if ( (min_vertex_(i) > b2.max_vertex_(i) + epsilon) ||
    			 (b2.min_vertex_(i)  > max_vertex_(i) + epsilon ) ) return false;
    	}
    	return true;
    }

    /**
     * Returns true if projection of the box to @p axis is an interval
     * less then (with tolerance) to given @p value.
     */
    bool projection_lt(unsigned int axis, double value) const
    {
    	return max_vertex_(axis) + epsilon < value;
    }

    /**
     * Returns true if projection of the box to @p axis is an interval
     * greater then (with tolerance) to given @p value.
     */
    bool projection_gt(unsigned int axis, double value) const
    {
    	return min_vertex_(axis) - epsilon > value;
    }

    /**
     * Split box into two boxes along @p axis by the
     * plane going through @p splitting_point on the axis.
     */
    void split(unsigned int axis, double splitting_point,
    		BoundingBox &left, BoundingBox &right ) const
    {
    	ASSERT_LESS(axis , dimension);
    	if (min_vertex_[axis] <= splitting_point && splitting_point <= max_vertex_[axis] ) {
    	   	left = *this;
    	   	right = *this;
    	 	left.max_vertex_[axis] = splitting_point;
    		right.min_vertex_[axis] = splitting_point;
    	} else {
    		THROW( ExcSplitting() << EI_interval_left(min_vertex_[axis])
    							  << EI_interval_right(max_vertex_[axis])
    							  << EI_split_point(splitting_point) );
    	}
    }

    /**
     * Expand bounding box to contain also given @p point.
     */
    void expand(const Point &point) {
		for(unsigned int j=0; j<dimension; j++) {
			min_vertex_(j) = std::min( min_vertex_[j], point[j] );
			max_vertex_(j) = std::max( max_vertex_[j], point[j] );
		}
    }

    /**
     * Return index of the axis in which the box has longest projection.
     */
   unsigned char longest_axis() const {
	   auto diff=max_vertex_ - min_vertex_;
	   return (diff[1] > diff[0])
			   	   ?  ( diff[2] > diff[1] ? 2 : 1 )
			   	   :  ( diff[2] > diff[0] ? 2 : 0 );
   }

private:
    /// stabilization parameter
    static const double epsilon;
    /// minimal coordinates of bounding box
    Point min_vertex_;
    /// maximal coordinates of bounding box
    Point max_vertex_;
};

/// Overloads output operator for box.
inline ostream &operator<<(ostream &stream, const BoundingBox &box) {
	stream << "Box("
		   << box.min(0) << " "
		   << box.min(1) << " "
		   << box.min(2) << "; "
		   << box.max(0) << " "
		   << box.max(1) << " "
		   << box.max(2) << ")";
	return stream;
}

#endif /* BOX_ELEMENT_HH_ */
