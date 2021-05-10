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
 *
 * TODO: remove virtual methods from Points classes, need changes also in Field::operator().
 */
class PointBase {
public:
    /// Default constructor
	PointBase()
    : local_point_idx_(0), elm_cache_map_(nullptr) {}

    /// Constructor
	PointBase(const ElementCacheMap *elm_cache_map, unsigned int loc_point_idx)
    : local_point_idx_(loc_point_idx), elm_cache_map_(elm_cache_map) {}

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

    /// Iterates to next point.
    void inc() {
    	this->local_point_idx_++;
    }

    /// Local coordinates within element
    //template<unsigned int dim>
    //inline arma::vec::fixed<dim> loc_coords() const {
    //    return this->eval_points()->local_point<dim>( this->eval_point_idx() );
    //}

    // Global coordinates within element
    //arma::vec3 coords() const;

protected:
    /// Index of the local point in the integral object.
    unsigned int local_point_idx_;
    /// Index of element in the patch.
    unsigned int elem_patch_idx_;
    /// Pointer ElementCacheMap needed for point evaluation.
    const ElementCacheMap* elm_cache_map_;
};


/**
 * @brief Point accessor allow iterate over bulk quadrature points defined in local element coordinates.
 */
class BulkPoint : public PointBase {
public:
    /// Default constructor
	BulkPoint()
    : PointBase() {}

    /// Constructor
	BulkPoint(const ElementCacheMap *elm_cache_map, PatchCacheLoc cache_pos)
    : PointBase(elm_cache_map, cache_pos.i_ep_) {
	    this->elem_patch_idx_ = cache_pos.i_elm_;
	}

    /// Comparison of accessors.
    bool operator==(const BulkPoint& other) {
        return (elem_patch_idx_ == other.elem_patch_idx_) && (local_point_idx_ == other.local_point_idx_);
    }

    /// Return index in EvalPoints object
    inline unsigned int eval_point_idx() const {
        return local_point_idx_;
    }
};


/**
 * @brief General point accessor allow iterate over quadrature points of given side defined in local element coordinates.
 *
 * Common ancestor of all side points classes (Edge-, Coupling-, BoundaryPoint)
 */
class SidePoint : public PointBase {
public:
    /// Default constructor
	SidePoint()
    : PointBase() {}

    /// Constructor
	SidePoint(DHCellSide cell_side, const ElementCacheMap *elm_cache_map, unsigned int local_point_idx)
    : PointBase(elm_cache_map, local_point_idx), side_idx_(cell_side.side_idx()),
	  permutation_idx_( cell_side.element()->permutation_idx( side_idx_ ) ) {
	    this->elem_patch_idx_ = this->elm_cache_map_->position_in_cache(cell_side.element().mesh_idx());
	}

    // Index of permutation
    inline unsigned int permutation_idx() const {
        return permutation_idx_;
    }

protected:
    /// Index of side in element
    unsigned int side_idx_;
    /// Permutation index corresponding with DHCellSide
    unsigned int permutation_idx_;
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
    EdgePoint(DHCellSide cell_side, const ElementCacheMap *elm_cache_map, const EdgeIntegral *edge_integral, unsigned int local_point_idx)
    : SidePoint(cell_side, elm_cache_map, local_point_idx), integral_(edge_integral) {}

    /// Return corresponds EdgePoint of neighbour side of same dimension (computing of side integrals).
    EdgePoint point_on(DHCellSide edg_side) const;

    /// Comparison of accessors.
    bool operator==(const EdgePoint& other) {
        return (elem_patch_idx_ == other.elem_patch_idx_) && (local_point_idx_ == other.local_point_idx_);
    }

    /// Return index in EvalPoints object
    inline unsigned int eval_point_idx() const;

private:
    /// Pointer to edge point integral
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
    CouplingPoint(DHCellSide cell_side, const ElementCacheMap *elm_cache_map, const CouplingIntegral *coupling_integral, unsigned int local_point_idx)
    : SidePoint(cell_side, elm_cache_map, local_point_idx), integral_(coupling_integral) {}

    /// Return corresponds EdgePoint of neighbour side of same dimension (computing of side integrals).
    BulkPoint lower_dim(DHCellAccessor cell_lower) const;

    /// Return index in EvalPoints object
    inline unsigned int eval_point_idx() const;

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
    BoundaryPoint(DHCellSide cell_side, const ElementCacheMap *elm_cache_map, const BoundaryIntegral *bdr_integral, unsigned int local_point_idx)
    : SidePoint(cell_side, elm_cache_map, local_point_idx), integral_(bdr_integral) {}

    /// Return corresponds BulkPoint on boundary element.
    BulkPoint point_bdr(ElementAccessor<3> bdr_elm) const;

    /// Return index in EvalPoints object
    inline unsigned int eval_point_idx() const;

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

    /// Constructor of bulk (n_permutations==0) or side subset
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
    /// Dimension of points
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
	BulkIntegral(std::shared_ptr<EvalPoints> eval_points, unsigned int dim)
	 : BaseIntegral(eval_points, dim), subset_index_(eval_points->n_subsets(dim)) {}

    /// Destructor
    ~BulkIntegral();

    /// Return index of data block according to subset in EvalPoints object
    inline int get_subset_idx() const {
        return subset_index_;
    }

    /// Returns range of bulk local points for appropriate cell accessor
    inline Range< BulkPoint > points(unsigned int element_patch_idx, const ElementCacheMap *elm_cache_map) const {
        auto bgn_it = make_iter<BulkPoint>( BulkPoint(elm_cache_map, PatchCacheLoc(element_patch_idx, eval_points_->subset_begin(dim_, subset_index_)) ));
        auto end_it = make_iter<BulkPoint>( BulkPoint(elm_cache_map, PatchCacheLoc(element_patch_idx, eval_points_->subset_end(dim_, subset_index_)) ));
        return Range<BulkPoint>(bgn_it, end_it);
    }

private:
    /// Index of data block according to subset in EvalPoints object.
    unsigned int subset_index_;
};

/**
 * Integral class of side points, allows assemblation of element - element fluxes.
 */
class EdgeIntegral : public BaseIntegral, public std::enable_shared_from_this<EdgeIntegral> {
public:
    /// Default constructor
	EdgeIntegral() : BaseIntegral(), perm_indices_(nullptr), n_permutations_(0) {}

    /// Constructor of edge integral
	EdgeIntegral(std::shared_ptr<EvalPoints> eval_points, unsigned int dim, unsigned int n_permutations, unsigned int points_per_side);

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

    /// Returns range of side local points for appropriate cell side accessor
    inline Range< EdgePoint > points(const DHCellSide &cell_side, const ElementCacheMap *elm_cache_map) const {
        unsigned int begin_idx = eval_points_->subset_begin(dim_, subset_index_);
        unsigned int end_idx = eval_points_->subset_end(dim_, subset_index_);
        unsigned int points_per_side = (end_idx - begin_idx) / this->n_sides();
        auto bgn_it = make_iter<EdgePoint>( EdgePoint(cell_side, elm_cache_map, this, 0 ) );
        auto end_it = make_iter<EdgePoint>( EdgePoint(cell_side, elm_cache_map, this, points_per_side ) );
        return Range<EdgePoint>(bgn_it, end_it);
    }

    /// Returns structure of permutation indices.
    inline int perm_idx_ptr(uint i_side, uint i_perm, uint i_point) const {
    	return perm_indices_[i_side][i_perm][i_point];
    }

private:
    /// Index of data block according to subset in EvalPoints object.
    unsigned int subset_index_;
    /// Indices to EvalPoints for different sides and permutations reflecting order of points.
    unsigned int*** perm_indices_;
    /// Number of sides (value 0 indicates bulk set)
    unsigned int n_sides_;
    /// Number of permutations (value 0 indicates bulk set)
    unsigned int n_permutations_;

    friend class EvalPoints;
    friend class EdgePoint;
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

    /// Returns range of side local points for appropriate cell side accessor
    inline Range< CouplingPoint > points(const DHCellSide &cell_side, const ElementCacheMap *elm_cache_map) const {
        unsigned int begin_idx = eval_points_->subset_begin(dim_, get_subset_high_idx());
        unsigned int end_idx = eval_points_->subset_end(dim_, get_subset_high_idx());
        unsigned int points_per_side = (end_idx - begin_idx) / edge_integral_->n_sides();
        auto bgn_it = make_iter<CouplingPoint>( CouplingPoint(cell_side, elm_cache_map, this, 0 ) );
        auto end_it = make_iter<CouplingPoint>( CouplingPoint(cell_side, elm_cache_map, this, points_per_side ) );
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

    /// Returns range of bulk local points for appropriate cell accessor
    inline Range< BoundaryPoint > points(const DHCellSide &cell_side, const ElementCacheMap *elm_cache_map) const {
        unsigned int begin_idx = eval_points_->subset_begin(dim_, get_subset_high_idx());
        unsigned int end_idx = eval_points_->subset_end(dim_, get_subset_high_idx());
        unsigned int points_per_side = (end_idx - begin_idx) / edge_integral_->n_sides();
        auto bgn_it = make_iter<BoundaryPoint>( BoundaryPoint(cell_side, elm_cache_map, this, 0 ) );
        auto end_it = make_iter<BoundaryPoint>( BoundaryPoint(cell_side, elm_cache_map, this, points_per_side ) );
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

inline unsigned int EdgePoint::eval_point_idx() const {
    return integral_->perm_idx_ptr(side_idx_, permutation_idx_, local_point_idx_);
}


inline unsigned int CouplingPoint::eval_point_idx() const {
    return integral_->edge_integral_->perm_idx_ptr(side_idx_, permutation_idx_, local_point_idx_);
}


inline unsigned int BoundaryPoint::eval_point_idx() const {
    return integral_->edge_integral_->perm_idx_ptr(side_idx_, permutation_idx_, local_point_idx_);
}

#endif /* EVAL_SUBSET_HH_ */
