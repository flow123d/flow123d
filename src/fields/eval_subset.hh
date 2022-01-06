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
 * @file    eval_subset.hh
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

#ifndef EVAL_SUBSET_HH_
#define EVAL_SUBSET_HH_

#include <memory>
#include <armadillo>
#include "fields/eval_points.hh"
#include "fields/field_value_cache.hh"
#include "mesh/range_wrapper.hh"
#include "mesh/accessors.hh"
#include "fem/dh_cell_accessor.hh"


class Side;
class BulkIntegral;
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

    /// Iterates to next point.
    void inc() {
    	this->local_point_idx_++;
    }

    /// Comparison of accessors.
    bool operator==(const BulkPoint& other) {
        ASSERT_EQ_DBG(elem_patch_idx_, other.elem_patch_idx_);
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

    /// Return index in EvalPoints object
    inline unsigned int eval_point_idx() const {
        return side_begin_ + local_point_idx_;
    }

    /// Comparison of accessors.
    bool operator==(const SidePoint& other) {
        ASSERT_EQ_DBG(elem_patch_idx_, other.elem_patch_idx_);
        ASSERT_EQ_DBG(side_begin_, other.side_begin_);
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
    inline EdgePoint(BulkPoint bulk, const EdgeIntegral *edge_integral, uint side_begin);

    /// Return corresponds EdgePoint of neighbour side of same dimension (computing of side integrals).
    inline EdgePoint point_on(const DHCellSide &edg_side) const;

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

    inline CouplingPoint(BulkPoint bulk, const CouplingIntegral *coupling_integral, uint side_begin);

    /// Return corresponds EdgePoint of neighbour side of same dimension (computing of side integrals).
    inline BulkPoint lower_dim(DHCellAccessor cell_lower) const;

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
    inline BoundaryPoint(BulkPoint bulk, const BoundaryIntegral *bdr_integral, uint side_begin);

    /// Return corresponds BulkPoint on boundary element.
    inline BulkPoint point_bdr(ElementAccessor<3> bdr_elm) const;

    /// Comparison of accessors.
    bool operator==(const BoundaryPoint& other) {
        return (elem_patch_idx_ == other.elem_patch_idx_) && (local_point_idx_ == other.local_point_idx_);
    }

private:
    /// Pointer to edge point set
    const BoundaryIntegral *integral_;
};



/**
 * Base integral class holds common data members and methods.
 */
class BaseIntegral {
public:
    /// Default constructor
	BaseIntegral() : eval_points_(nullptr), dim_(0) {}

    /// Constructor of bulk or side subset
	BaseIntegral(std::shared_ptr<EvalPoints> eval_points, unsigned int dim)
	 : eval_points_(eval_points), dim_(dim) {}

    /// Destructor
    virtual ~BaseIntegral();

    /// Getter of eval_points
    std::shared_ptr<EvalPoints> eval_points() const {
        return eval_points_;
    }

    /// Returns dimension.
    unsigned int dim() const {
    	return dim_;
    }
protected:
    /// Pointer to EvalPoints
    std::shared_ptr<EvalPoints> eval_points_;
    /// Dimension of the cell on which points are placed
    unsigned int dim_;
};

/**
 * Integral class of bulk points, allows assemblation of volume integrals.
 */
class BulkIntegral : public BaseIntegral, public std::enable_shared_from_this<BulkIntegral> {
public:
    /// Default constructor
	BulkIntegral() : BaseIntegral() {}

    /// Constructor of bulk integral
	BulkIntegral(std::shared_ptr<EvalPoints> eval_points, unsigned int dim, uint i_subset)
	 : BaseIntegral(eval_points, dim), subset_index_(i_subset)
	{
	    begin_idx_ = eval_points_->subset_begin(dim_, subset_index_);
	    end_idx_ = eval_points_->subset_end(dim_, subset_index_);
	}

    /// Destructor
    ~BulkIntegral();

    /// Return index of data block according to subset in EvalPoints object
    inline int get_subset_idx() const {
        return subset_index_;
    }


    /// Returns range of bulk local points for appropriate cell accessor
    inline Range< BulkPoint > points(unsigned int element_patch_idx, const ElementCacheMap *elm_cache_map) const {
        auto bgn_it = make_iter<BulkPoint>( BulkPoint(elm_cache_map, element_patch_idx, begin_idx_));
        auto end_it = make_iter<BulkPoint>( BulkPoint(elm_cache_map, element_patch_idx, end_idx_));
        return Range<BulkPoint>(bgn_it, end_it);
    }

protected:
    /// Index of data block according to subset in EvalPoints object.
    unsigned int subset_index_;
    uint begin_idx_;
    uint end_idx_;

};

/**
 * Integral class of side points, allows assemblation of element - element fluxes.
 */
class EdgeIntegral : public BaseIntegral, public std::enable_shared_from_this<EdgeIntegral> {
public:
    /// Default constructor
	EdgeIntegral() : BaseIntegral()
    {
	    ASSERT(false);
    }

    /// Constructor of edge integral
	EdgeIntegral(std::shared_ptr<EvalPoints> eval_points, unsigned int dim, uint i_subset);

    /// Destructor
    ~EdgeIntegral();

    /// Getter of n_sides
    inline unsigned int n_sides() const {
        return n_sides_;
    }

    /// Return index of data block according to subset in EvalPoints object
    inline int get_subset_idx() const {
        return subset_index_;
    }

    inline uint side_begin(const DHCellSide &cell_side) const {
        return begin_idx_ + cell_side.side_idx() * n_points_per_side_;
    }

    /// Returns range of side local points for appropriate cell side accessor
    inline Range< EdgePoint > points(const DHCellSide &cell_side, const ElementCacheMap *elm_cache_map) const {
        ASSERT_EQ_DBG(cell_side.dim(), dim_);
        //DebugOut() << "points per side: " << n_points_per_side_;
        uint element_patch_idx = elm_cache_map->position_in_cache(cell_side.element().idx());
        uint begin_idx = side_begin(cell_side);
        auto bgn_it = make_iter<EdgePoint>( EdgePoint(
                BulkPoint(elm_cache_map, element_patch_idx, 0), this, begin_idx));
        auto end_it = make_iter<EdgePoint>( EdgePoint(
                BulkPoint(elm_cache_map, element_patch_idx, n_points_per_side_), this, begin_idx));
        return Range<EdgePoint>(bgn_it, end_it);
    }


private:
    unsigned int subset_index_;
    uint begin_idx_;

    /// Number of sides (value 0 indicates bulk set)
    unsigned int n_sides_;
    /// Number of points. TODO: pass this to the constructor, avoid extraction from the eval_points
    uint n_points_per_side_;

    friend class EvalPoints;
    friend class EdgePoint;
    friend class CouplingIntegral;
    friend class BoundaryIntegral;
};

/**
 * Integral class of neighbour points, allows assemblation of element - side fluxes.
 *
 * Dimension corresponds with element of higher dim.
 */
class CouplingIntegral : public BaseIntegral, public std::enable_shared_from_this<CouplingIntegral> {
public:
    /// Default constructor
	CouplingIntegral() : BaseIntegral() {}

    /// Constructor of ngh integral
	CouplingIntegral(std::shared_ptr<EdgeIntegral> edge_integral, std::shared_ptr<BulkIntegral> bulk_integral);

    /// Destructor
    ~CouplingIntegral();

    /// Return index of data block according to subset of higher dim in EvalPoints object
    inline int get_subset_high_idx() const {
        return edge_integral_->get_subset_idx();
    }

    /// Return index of data block according to subset of lower dim in EvalPoints object
    inline int get_subset_low_idx() const {
        return bulk_integral_->get_subset_idx();
    }

    inline uint bulk_begin() const {
        return eval_points_->subset_begin(dim_-1, bulk_integral_->get_subset_idx());
    }

    /// Returns range of side local points for appropriate cell side accessor
    inline Range< CouplingPoint > points(const DHCellSide &cell_side, const ElementCacheMap *elm_cache_map) const {
        ASSERT_EQ_DBG(cell_side.dim(), dim_);
        uint element_patch_idx = elm_cache_map->position_in_cache(cell_side.element().idx());
        uint begin_idx = edge_integral_->side_begin(cell_side);
        auto bgn_it = make_iter<CouplingPoint>( CouplingPoint(
                BulkPoint(elm_cache_map, element_patch_idx, 0), this, begin_idx) );
        auto end_it = make_iter<CouplingPoint>( CouplingPoint(
                BulkPoint(elm_cache_map, element_patch_idx, edge_integral_->n_points_per_side_), this, begin_idx) );;
        return Range<CouplingPoint>(bgn_it, end_it);
    }

private:
    /// Integral according to side subset part (element of higher dim) in EvalPoints object.
    std::shared_ptr<EdgeIntegral> edge_integral_;
    /// Integral according to bulk subset part (element of lower dim) in EvalPoints object.
    std::shared_ptr<BulkIntegral> bulk_integral_;

    friend class CouplingPoint;
};

/**
 * Integral class of boundary points, allows assemblation of fluxes between sides and neighbouring boundary elements.
 */
class BoundaryIntegral : public BaseIntegral, public std::enable_shared_from_this<BoundaryIntegral> {
public:
    /// Default constructor
    BoundaryIntegral() : BaseIntegral() {}

    /// Constructor of bulk subset
    BoundaryIntegral(std::shared_ptr<EdgeIntegral> edge_integral, std::shared_ptr<BulkIntegral> bulk_integral);

    /// Destructor
    ~BoundaryIntegral();

    /// Return index of data block according to subset of higher dim in EvalPoints object
    inline int get_subset_high_idx() const {
        return edge_integral_->get_subset_idx();
    }

    /// Return index of data block according to subset of lower dim (boundary) in EvalPoints object
    inline int get_subset_low_idx() const {
        return bulk_integral_->get_subset_idx();
    }

    inline uint bulk_begin() const {
      //  DebugOut().fmt("edge_begin: {} bdr_begin: {}",
      //          eval_points_->subset_begin(dim_, edge_integral_->get_subset_idx()),
      //          eval_points_->subset_begin(dim_-1, bulk_integral_->get_subset_idx()));
        return eval_points_->subset_begin(dim_-1, bulk_integral_->get_subset_idx());
    }

    /// Returns range of bulk local points for appropriate cell accessor
    inline Range< BoundaryPoint > points(const DHCellSide &cell_side, const ElementCacheMap *elm_cache_map) const {
        ASSERT_EQ_DBG(cell_side.dim(), dim_);
        uint element_patch_idx = elm_cache_map->position_in_cache(cell_side.element().idx());
        uint begin_idx = edge_integral_->side_begin(cell_side);
        auto bgn_it = make_iter<BoundaryPoint>( BoundaryPoint(
                BulkPoint(elm_cache_map, element_patch_idx, 0), this, begin_idx) );
        auto end_it = make_iter<BoundaryPoint>( BoundaryPoint(
                BulkPoint(elm_cache_map, element_patch_idx, edge_integral_->n_points_per_side_), this, begin_idx) );;
        return Range<BoundaryPoint>(bgn_it, end_it);
    }

private:
    /// Integral according to higher dim (bulk) element subset part in EvalPoints object.
    std::shared_ptr<EdgeIntegral> edge_integral_;
    /// Integral according to kower dim (boundary) element subset part in EvalPoints object.
    std::shared_ptr<BulkIntegral> bulk_integral_;

    friend class BoundaryPoint;
};


/******************************************************************************
 * Implementation of inlined methods
 */

/// Constructor
EdgePoint::EdgePoint(BulkPoint bulk, const EdgeIntegral *edge_integral, uint side_begin)
: SidePoint(bulk, side_begin),
  integral_(edge_integral)
{}

inline EdgePoint EdgePoint::point_on(const DHCellSide &edg_side) const {
    uint element_patch_idx = elm_cache_map_->position_in_cache(edg_side.element().idx());
    uint side_begin = integral_->side_begin(edg_side);
    return EdgePoint(BulkPoint(elm_cache_map_, element_patch_idx, local_point_idx_),
            integral_, side_begin);
}

//******************************************************************************
CouplingPoint::CouplingPoint(BulkPoint bulk, const CouplingIntegral *coupling_integral, uint side_begin)
: SidePoint(bulk, side_begin),
  integral_(coupling_integral)
{}


inline BulkPoint CouplingPoint::lower_dim(DHCellAccessor cell_lower) const {
    unsigned int i_elm = elm_cache_map_->position_in_cache(cell_lower.elm().idx());
    unsigned int i_ep = integral_->bulk_begin() + local_point_idx_;
    return BulkPoint(elm_cache_map_, i_elm, i_ep);
}



//******************************************************************************
BoundaryPoint::BoundaryPoint(BulkPoint bulk, const BoundaryIntegral *bdr_integral, uint side_begin)
: SidePoint(bulk, side_begin),
  integral_(bdr_integral)
{}



inline BulkPoint BoundaryPoint::point_bdr(ElementAccessor<3> bdr_elm) const {
    unsigned int i_elm = elm_cache_map_->position_in_cache(bdr_elm.idx(), true);
    unsigned int i_ep = integral_->bulk_begin() + local_point_idx_;
    //DebugOut() << "begin:" << integral_->bulk_begin() << "iloc " << local_point_idx_;
    return BulkPoint(elm_cache_map_, i_elm, i_ep);
}

#endif /* EVAL_SUBSET_HH_ */
