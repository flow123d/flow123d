/*!
 *
﻿ * Copyright (C) 2015 Technical University of Liberec.  All rights reserved.
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
 * @file    integral_points.hh
 * @brief
 * @author  David Flanderka
 *
 * TODO (readability, optimization):
 * - EdgePoint::point_on is too general, can fail for non-matching sides,
 *   depends on a map;
 *   Iterator over edge sides should know index of side on the edge or
 *   we should directly iterate over pairs of sides and iterate points over both sides
 *   without mapping.
 * - similarly for Coupling and boundary points
 * - Points should just hold necessary indices, without reference to complex classes,
 *   these points only can be used as indices to fields to get appropriate value in the field cache
 */

#ifndef INTEGRAL_POINTS_HH_
#define INTEGRAL_POINTS_HH_

#include <memory>
#include <armadillo>
#include "fem/eval_points.hh"
#include "fem/element_cache_map.hh"
#include "mesh/accessors.hh"
#include "fem/dh_cell_accessor.hh"


class Side;
class EdgeIntegral;
class CouplingIntegral;
class BoundaryIntegral;


/**
 * @brief Base point accessor class.
 */
/**
 * @brief Point accessor allow iterate over bulk quadrature points defined in local element coordinates.
 */

class BulkPoint {
public:
    /// Default constructor
    BulkPoint()
    {}

    /// Constructor
    BulkPoint(const ElementCacheMap *elm_cache_map, uint elem_idx, uint loc_point_idx)
    : elm_cache_map_(elm_cache_map), elem_patch_idx_(elem_idx), local_point_idx_(loc_point_idx)
	{}

    /// Getter of EvalPoints object.
    inline std::shared_ptr<EvalPoints> eval_points() const {
        ASSERT_PTR(elm_cache_map_).error("Invalid point.\n");
        return elm_cache_map_->eval_points();
    }

    // Getter of ElementCacheMap object.
    inline const ElementCacheMap *elm_cache_map() const {
        return elm_cache_map_;
    }


	// Getter of element patch index.
    inline unsigned int elem_patch_idx() const {
        return elem_patch_idx_;
    }

    /// Return index in EvalPoints object
    inline unsigned int eval_point_idx() const {
        return local_point_idx_;
    }

    /// Return index in ElementCacheMap
    inline unsigned int value_cache_idx() const {
        return elm_cache_map_->element_eval_point(elem_patch_idx_, local_point_idx_);
    }

    /// Iterates to next point.
    void inc() {
    	this->local_point_idx_++;
    }

    /// Comparison of accessors.
    bool operator==(const BulkPoint& other) {
        ASSERT_EQ(elem_patch_idx_, other.elem_patch_idx_);
        return (local_point_idx_ == other.local_point_idx_);
    }


protected:
    const ElementCacheMap *elm_cache_map_;
    /// Index of element in the patch.
    unsigned int elem_patch_idx_;
    /// Index of the local point in the integral object.
    unsigned int local_point_idx_;
};




/**
 * @brief General point a+ side_begin_ + ccessor allow iterate over quadrature points of given side defined in local element coordinates.
 *
 * Common ancestor of all side points classes (Edge-, Coupling-, BoundaryPoint)
 */
class SidePoint : public BulkPoint {
public:
    /// Default constructor
	SidePoint()
    : BulkPoint() {}

	/// Constructor
	SidePoint(BulkPoint bulk, uint side_begin)
	: BulkPoint(bulk), side_begin_(side_begin)
	{
	    //DebugOut().fmt("begin: {} sidx: {}", side_begin_, local_point_idx_);
	}


    /// Constructor
	inline SidePoint(DHCellSide cell_side, const ElementCacheMap *elm_cache_map,
	        const EdgeIntegral *edge_integral, unsigned int local_point_idx);

    /// Return local index in quadrature. Temporary method - intermediate step in implementation of PatcFEValues.
    inline unsigned int local_point_idx() const {
        return local_point_idx_;
    }

    /// Return index in EvalPoints object
    inline unsigned int eval_point_idx() const {
        return side_begin_ + local_point_idx_;
    }

    /// Return index in ElementCacheMap
    inline unsigned int value_cache_idx() const {
        return elm_cache_map_->element_eval_point(elem_patch_idx_, eval_point_idx());
    }

    /// Comparison of accessors.
    bool operator==(const SidePoint& other) {
        ASSERT_EQ(elem_patch_idx_, other.elem_patch_idx_);
        ASSERT_EQ(side_begin_, other.side_begin_);
        return (local_point_idx_ == other.local_point_idx_);
    }


protected:
    //// local_point_idx_ here have meaning of index within the side.

    /// Index of side in element
    unsigned int side_begin_;

};


/**
 * @brief Point accessor allow iterate over quadrature points of given side defined in local element coordinates.
 */
class EdgePoint : public SidePoint {
public:
    /// Default constructor
    EdgePoint()
    : SidePoint() {}

    /// Constructor
    inline EdgePoint(BulkPoint bulk, const EdgeIntegral *edge_integral, uint side_begin)
    : SidePoint(bulk, side_begin),
      integral_(edge_integral)
    {}

    /// Return corresponds EdgePoint of neighbour side of same dimension (computing of side integrals).
    EdgePoint point_on(const DHCellSide &edg_side) const;

    /// Comparison of accessors.
    bool operator==(const EdgePoint& other) {
        return (elem_patch_idx_ == other.elem_patch_idx_) && (local_point_idx_ == other.local_point_idx_);
    }
private:
    const EdgeIntegral *integral_;
};


/**
 * @brief Point accessor allow iterate over quadrature points of given side defined in local element coordinates.
 */
class CouplingPoint : public SidePoint {
public:
    /// Default constructor
    CouplingPoint()
    : SidePoint() {}

    /// Constructor
    inline CouplingPoint(BulkPoint bulk, const CouplingIntegral *coupling_integral, uint side_begin)
    : SidePoint(bulk, side_begin),
      integral_(coupling_integral)
    {}

    /// Return corresponds EdgePoint of neighbour side of same dimension (computing of side integrals).
    BulkPoint lower_dim(DHCellAccessor cell_lower) const;

    /// Comparison of accessors.
    bool operator==(const CouplingPoint& other) {
        return (elem_patch_idx_ == other.elem_patch_idx_) && (local_point_idx_ == other.local_point_idx_);
    }

private:
    /// Pointer to edge point set
    const CouplingIntegral *integral_;
};

/**
 * @brief Point accessor allow iterate over quadrature points of given side defined in local element coordinates.
 */
class BoundaryPoint : public SidePoint {
public:
    /// Default constructor
    BoundaryPoint()
    : SidePoint() {}

    /// Constructor
    inline BoundaryPoint(BulkPoint bulk, const BoundaryIntegral *bdr_integral, uint side_begin)
    : SidePoint(bulk, side_begin),
      integral_(bdr_integral)
    {}

    /// Return corresponds BulkPoint on boundary element.
    BulkPoint point_bdr(ElementAccessor<3> bdr_elm) const;

    /// Comparison of accessors.
    bool operator==(const BoundaryPoint& other) {
        return (elem_patch_idx_ == other.elem_patch_idx_) && (local_point_idx_ == other.local_point_idx_);
    }

private:
    /// Pointer to edge point set
    const BoundaryIntegral *integral_;
};


#endif /* INTEGRAL_POINTS_HH_ */
