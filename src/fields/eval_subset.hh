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
#include "fem/dh_cell_accessor.hh"


class Side;
class BulkPoint;
class EdgePoint;


/**
 * Base integral class holds common data members and methods.
 */
class BaseIntegral {
public:
	TYPEDEF_ERR_INFO(EI_ElementIdx, unsigned int);
	DECLARE_EXCEPTION(ExcElementNotInCache,
	        << "Element of Idx: " << EI_ElementIdx::val << " is not stored in 'Field value data cache'.\n"
			   << "Value can't be computed.\n");

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
    int get_subset_idx() const {
        return subset_index_;
    }

    /// Returns range of bulk local points for appropriate cell accessor
    Range< BulkPoint > points(const DHCellAccessor &cell, const ElementCacheMap *elm_cache_map) const;

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
    unsigned int n_sides() const {
        return n_sides_;
    }

    /// Return index of data block according to subset in EvalPoints object
    int get_subset_idx() const {
        return subset_index_;
    }

    /// Returns range of side local points for appropriate cell side accessor
    Range< EdgePoint > points(const DHCellSide &cell_side, const ElementCacheMap *elm_cache_map) const;

    /// Returns structure of permutation indices.
    int perm_idx_ptr(uint i_side, uint i_perm, uint i_point) const {
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
    int get_subset_high_idx() const {
        return edge_integral_->get_subset_idx();
    }

    /// Return index of data block according to subset of lower dim in EvalPoints object
    int get_subset_low_idx() const {
        return bulk_integral_->get_subset_idx();
    }

    /// Returns range of bulk local points for appropriate cell accessor
    Range< BulkPoint > points(const DHCellAccessor &cell, const ElementCacheMap *elm_cache_map) const;

    /// Returns range of side local points for appropriate cell side accessor
    Range< EdgePoint > points(const DHCellSide &cell_side, const ElementCacheMap *elm_cache_map) const;

private:
    /// Integral according to side subset part (element of higher dim) in EvalPoints object.
    std::shared_ptr<EdgeIntegral> edge_integral_;
    /// Integral according to bulk subset part (element of lower dim) in EvalPoints object.
    std::shared_ptr<BulkIntegral> bulk_integral_;
};

/**
 * Integral class of boundary points, allows assemblation of fluxes between sides and neighbouring boundary elements.
 */
class BoundaryIntegral : public BaseIntegral {
public:
    /// Default constructor
    BoundaryIntegral() : BaseIntegral() {}

    /// Constructor of bulk subset
    BoundaryIntegral(std::shared_ptr<EdgeIntegral> edge_integral);

    /// Destructor
    ~BoundaryIntegral();

    /// Return index of data block according to subset in EvalPoints object
    int get_subset_idx() const {
        return edge_integral_->get_subset_idx();
    }

    /// Returns range of bulk local points for appropriate cell accessor
    Range< EdgePoint > points(const DHCellSide &cell_side, const ElementCacheMap *elm_cache_map) const;

private:
    /// Boundary integral according to edge integral (? but need own special data members and methods ?).
    std::shared_ptr<EdgeIntegral> edge_integral_;
};


/**
 * @brief Point accessor allow iterate over bulk quadrature points defined in local element coordinates.
 */
class BulkPoint {
public:
    /// Default constructor
	BulkPoint()
    : local_point_idx_(0) {}

    /// Constructor
	BulkPoint(DHCellAccessor dh_cell, const ElementCacheMap *elm_cache_map, std::shared_ptr<const BulkIntegral> bulk_integral, unsigned int loc_point_idx)
    : dh_cell_(dh_cell), integral_(bulk_integral), local_point_idx_(loc_point_idx), elm_cache_map_(elm_cache_map) {}

    /// Getter of BulkIntegral
    std::shared_ptr<const BulkIntegral> integral() const {
        return integral_;
    }

    /// Getter of EvalPoints
    std::shared_ptr<EvalPoints> eval_points() const {
        return integral_->eval_points();
    }

    /// Local coordinates within element
    template<unsigned int dim>
    inline arma::vec::fixed<dim> loc_coords() const {
        return this->eval_points()->local_point<dim>( local_point_idx_ );
    }

    // Global coordinates within element
    //arma::vec3 coords() const;

    /// Return index of element in data cache.
    unsigned int element_cache_index() const {
        return dh_cell_.element_cache_index();
    }

    /// Return DH cell accessor.
    DHCellAccessor dh_cell() const {
        return dh_cell_;
    }

    // Index of permutation
    inline const ElementCacheMap *elm_cache_map() const {
        return elm_cache_map_;
    }

    /// Return index in EvalPoints object
    unsigned int eval_point_idx() const {
        return local_point_idx_;
    }

    /// Iterates to next point.
    void inc() {
    	local_point_idx_++;
    }

    /// Comparison of accessors.
    bool operator==(const BulkPoint& other) {
    	return (dh_cell_ == other.dh_cell_) && (local_point_idx_ == other.local_point_idx_);
    }

private:
    /// DOF handler accessor of element.
    DHCellAccessor dh_cell_;
    /// Pointer to bulk integral.
    std::shared_ptr<const BulkIntegral> integral_;
    /// Index of the local point in bulk point set.
    unsigned int local_point_idx_;
    /// Pointer ElementCacheMap needed for point evaluation.
    const ElementCacheMap* elm_cache_map_;
};


/**
 * @brief Point accessor allow iterate over quadrature points of given side defined in local element coordinates.
 */
class EdgePoint {
public:
    /// Default constructor
	EdgePoint()
    : local_point_idx_(0), elm_cache_map_(nullptr) {}

    /// Constructor
	EdgePoint(DHCellSide cell_side, const ElementCacheMap *elm_cache_map, std::shared_ptr<const EdgeIntegral> edge_integral, unsigned int local_point_idx)
    : cell_side_(cell_side), integral_(edge_integral), local_point_idx_(local_point_idx),
	  permutation_idx_( cell_side.element()->permutation_idx( cell_side_.side_idx() ) ), elm_cache_map_(elm_cache_map) {}

    /// Getter of EdgeIntegral
    std::shared_ptr<const EdgeIntegral> integral() const {
        return integral_;
    }

    /// Getter of evaluation points
    std::shared_ptr<EvalPoints> eval_points() const {
        return integral_->eval_points();
    }

    // Local coordinates within element
    template<unsigned int dim>
    inline arma::vec::fixed<dim> loc_coords() const {
        return this->eval_points()->local_point<dim>( this->eval_point_idx() );
    }

    // Global coordinates within element
    //arma::vec3 coords() const;

    /// Return index of element in data cache.
    unsigned int element_cache_index() const {
        return cell_side_.cell().element_cache_index();
    }

    /// Return DH cell accessor.
    DHCellSide dh_cell_side() const {
        return cell_side_;
    }

    // Index of permutation
    unsigned int permutation_idx() const {
        return permutation_idx_;
    }

    // Index of permutation
    inline const ElementCacheMap *elm_cache_map() const {
        return elm_cache_map_;
    }

    /// Return index in EvalPoints object
    unsigned int eval_point_idx() const {
        return integral_->perm_idx_ptr(cell_side_.side_idx(), permutation_idx_, local_point_idx_);
    }

    /// Return corresponds EdgePoint of neighbour side of same dimension (computing of side integrals).
    EdgePoint permute(DHCellSide edg_side) const;

    /// Iterates to next point.
    void inc() {
    	local_point_idx_++;
    }

    /// Comparison of accessors.
    bool operator==(const EdgePoint& other) {
    	return (cell_side_ == other.cell_side_) && (local_point_idx_ == other.local_point_idx_);
    }

private:
    /// DOF handler accessor of element side.
    DHCellSide cell_side_;
    /// Pointer to edge point set
    std::shared_ptr<const EdgeIntegral> integral_;
    /// Index of the local point in the composed quadrature.
    unsigned int local_point_idx_;
    /// Permutation index corresponding with DHCellSide
    unsigned int permutation_idx_;
    /// Pointer ElementCacheMap needed for point evaluation.
    const ElementCacheMap* elm_cache_map_;
};



#endif /* EVAL_SUBSET_HH_ */
