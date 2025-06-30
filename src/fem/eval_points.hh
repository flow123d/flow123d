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
 * @file    eval_points.hh
 * @brief
 * @author  David Flanderka
 */

#ifndef EVAL_POINTS_HH_
#define EVAL_POINTS_HH_


#include <vector>
#include <memory>
#include <armadillo>
#include "fem/dh_cell_accessor.hh"
#include "tools/revertable_list.hh"
#include "mesh/range_wrapper.hh"
#include "quadrature/quadrature.hh"
#include "system/asserts.hh"
#include "system/armor.hh"

class Side;
class BulkIntegral;
class EdgeIntegral;
class CouplingIntegral;
class BoundaryIntegral;
template <int spacedim> class ElementAccessor;


/**
 * Helper structure holds data of cell (bulk) integral
 *
 * Data is specified by cell and subset index in EvalPoint object
 */
struct BulkIntegralData {
	/// Default constructor
    BulkIntegralData() {}

    /// Constructor with data mebers initialization
    BulkIntegralData(DHCellAccessor dhcell, unsigned int subset_idx)
    : cell(dhcell), subset_index(subset_idx) {}

    /// Copy constructor
    BulkIntegralData(const BulkIntegralData &other)
    : cell(other.cell), subset_index(other.subset_index) {}

    DHCellAccessor cell;          ///< Specified cell (element)
    unsigned int subset_index;    ///< Index (order) of subset in EvalPoints object
};

/**
 * Helper structure holds data of edge integral
 *
 * Data is specified by side and subset index in EvalPoint object
 */
struct EdgeIntegralData {
	/// Default constructor
	EdgeIntegralData()
	: edge_side_range(make_iter<DHEdgeSide, DHCellSide>( DHEdgeSide() ), make_iter<DHEdgeSide, DHCellSide>( DHEdgeSide() )) {}

    /// Copy constructor
	EdgeIntegralData(const EdgeIntegralData &other)
    : edge_side_range(other.edge_side_range), subset_index(other.subset_index) {}

    /// Constructor with data mebers initialization
	EdgeIntegralData(RangeConvert<DHEdgeSide, DHCellSide> range, unsigned int subset_idx)
    : edge_side_range(range), subset_index(subset_idx) {}

	RangeConvert<DHEdgeSide, DHCellSide> edge_side_range;   ///< Specified cell side (element)
    unsigned int subset_index;                              ///< Index (order) of subset in EvalPoints object
};

/**
 * Helper structure holds data of neighbour (coupling) integral
 *
 * Data is specified by cell, side and their subset indices in EvalPoint object
 */
struct CouplingIntegralData {
	/// Default constructor
   	CouplingIntegralData() {}

    /// Constructor with data mebers initialization
   	CouplingIntegralData(DHCellAccessor dhcell, unsigned int bulk_idx, DHCellSide dhside, unsigned int side_idx)
    : cell(dhcell), bulk_subset_index(bulk_idx), side(dhside), side_subset_index(side_idx) {}

    /// Copy constructor
   	CouplingIntegralData(const CouplingIntegralData &other)
    : cell(other.cell), bulk_subset_index(other.bulk_subset_index), side(other.side), side_subset_index(other.side_subset_index) {}

    DHCellAccessor cell;
    unsigned int bulk_subset_index;    ///< Index (order) of lower dim subset in EvalPoints object
    DHCellSide side;                   ///< Specified cell side (higher dim element)
    unsigned int side_subset_index;    ///< Index (order) of higher dim subset in EvalPoints object
};

/**
 * Helper structure holds data of boundary integral
 *
 * Data is specified by side and subset indices of side and appropriate boundary element in EvalPoint object
 */
struct BoundaryIntegralData {
	/// Default constructor
	BoundaryIntegralData() {}

    /// Constructor with data mebers initialization
	BoundaryIntegralData(unsigned int bdr_idx, DHCellSide dhside, unsigned int side_idx)
    : bdr_subset_index(bdr_idx), side(dhside), side_subset_index(side_idx) {}

    /// Copy constructor
	BoundaryIntegralData(const BoundaryIntegralData &other)
    : bdr_subset_index(other.bdr_subset_index), side(other.side), side_subset_index(other.side_subset_index) {}

	// We don't need hold ElementAccessor of boundary element, side.cond().element_accessor() provides it.
    unsigned int bdr_subset_index;     ///< Index (order) of subset on boundary element in EvalPoints object
    DHCellSide side;                   ///< Specified cell side (bulk element)
    unsigned int side_subset_index;    ///< Index (order) of subset on side of bulk element in EvalPoints object
};


/// Set of integral data of given dimension used in assemblation
struct IntegralData {
public:
    RevertableList<BulkIntegralData>       bulk_;      ///< Holds data for computing bulk integrals.
    RevertableList<EdgeIntegralData>       edge_;      ///< Holds data for computing edge integrals.
    RevertableList<CouplingIntegralData>   coupling_;  ///< Holds data for computing couplings integrals.
    RevertableList<BoundaryIntegralData>   boundary_;  ///< Holds data for computing boundary integrals.
};



template <class Integral>
class IntegralHashMap {
private:
    /// Return unique hash
    inline uint quad_hash(Quadrature *quad, unsigned int dim) const {
        return dim + 4*quad->size();
    }

    std::unordered_map< uint, std::shared_ptr<Integral> > integral_map_;

    /// Value-only iterator
    class ValueIterator {
        using InnerIterator = typename std::unordered_map< uint, std::shared_ptr<Integral> >::const_iterator;
        InnerIterator it;

    public:
        explicit ValueIterator(InnerIterator iter) : it(iter) {}

        const std::shared_ptr<Integral>& operator*() const { return it->second; }
        ValueIterator& operator++() { ++it; return *this; }
        bool operator!=(const ValueIterator& other) const { return it != other.it; }
    };
public:
    /// Constructor
	IntegralHashMap() {}

    /// Value-only iterable wrapper
    class ValueIterable {
        const std::unordered_map< uint, std::shared_ptr<Integral> >& integral_map_ref_;
    public:
        ValueIterable(const std::unordered_map< uint, std::shared_ptr<Integral> >& ref) : integral_map_ref_(ref) {}

        auto begin() const { return ValueIterator(integral_map_ref_.begin()); }
        auto end() const { return ValueIterator(integral_map_ref_.end()); }
    };

    /// Operator returning map values iterable
    ValueIterable operator()() const {
        return ValueIterable(integral_map_);
    }

    /// Checks if object of same quadrature exists, if not creates it and return it.
    std::shared_ptr<Integral> create_and_return(Quadrature *quad, unsigned int dim) {
    	uint hash = quad_hash(quad, dim);
        auto it = integral_map_.find(hash);

        if (it == integral_map_.end()) {
        	integral_map_[hash] = std::make_shared<Integral>(quad, dim);
            it = integral_map_.find(hash);
        }
        return it->second;
    }

    //// Define begin() and end() methods to allow iteration
    auto begin() { return integral_map_.begin(); }
    auto end() { return integral_map_.end(); }

    //// Constant modification of previous methods
    auto begin() const { return integral_map_.begin(); }
    auto end() const { return integral_map_.end(); }

    /// Return size of integral map
    unsigned int size() const {
    	return integral_map_.size();
    }

    /// Copy all items of other IntegralHashMap to this.
    void copy_items(IntegralHashMap<Integral> &other) {
        for (auto it=other.begin(); it!=other.end(); ++it) {
            auto own_it = integral_map_.find(it->first);
            ASSERT(own_it == integral_map_.end());
            integral_map_[it->first] = it->second;
        }
    }

};


/// Set of integral of given dimension necessary in assemblation
struct DimIntegrals {
	IntegralHashMap<BulkIntegral> bulk_;            ///< Bulk integrals of elements
	IntegralHashMap<EdgeIntegral> edge_;            ///< Edge integrals between elements of same dimensions
	IntegralHashMap<CouplingIntegral> coupling_;    ///< Coupling integrals between elements of dimensions dim and dim-1
	IntegralHashMap<BoundaryIntegral> boundary_;    ///< Boundary integrals betwwen side and boundary element of dim-1
};


/**
 * @brief Class holds local coordinations of evaluating points (bulk and sides)
 * specified by element dimension.
 */
class EvalPoints : public std::enable_shared_from_this<EvalPoints> {
public:
	/// Undefined dimension of new (empty) object
	static const unsigned int undefined_dim;

	/// Maximal number of hold subsets.
	static constexpr unsigned int max_subsets = 10;

	/// Maximal average number of points hold in subset.
	static const unsigned int max_subset_points = 30;

    /// Constructor
	EvalPoints();

    /// Return size of evaluation points object (number of points).
    inline unsigned int size(unsigned int dim) const {
        return dim_eval_points_[dim].size();
    }

    /// Return local coordinates of given local point and appropriate dim.
    template<unsigned int dim>
    inline arma::vec::fixed<dim> local_point(unsigned int local_point_idx) const {
    	ASSERT_GT(dim, 0).error("Dimension 0 not supported!\n");
        return dim_eval_points_[dim].local_point<dim>(local_point_idx);
    }

    /// Return begin index of appropriate subset data.
    inline int subset_begin(unsigned int dim, unsigned int idx) const {
        return dim_eval_points_[dim].subset_begin(idx);
    }

    /// Return end index of appropriate subset data.
    inline int subset_end(unsigned int dim, unsigned int idx) const {
        return dim_eval_points_[dim].subset_end(idx);
    }

    /// Return number of local points corresponding to subset.
    inline int subset_size(unsigned int dim, unsigned int idx) const {
        return dim_eval_points_[dim].subset_size(idx);
    }

    /// Return number of subsets.
    inline unsigned int n_subsets(unsigned int dim) const {
        return dim_eval_points_[dim].n_subsets();
    }

    /**
     * Registers point set from quadrature.
     * Returns an object referencing to the EvalPoints and list of its points.
     */
    template <unsigned int dim>
    unsigned int add_bulk(const Quadrature &);

    /// The same as add_bulk but for edge points on sides.
    template <unsigned int dim>
    unsigned int add_edge(const Quadrature &);

//    /// The same as add_bulk but for points between side points of element of dim and bulk points of element of dim-1.
//    template <unsigned int dim>
//    std::shared_ptr<CouplingIntegral> add_coupling(const Quadrature &);
//
//    /// The same as add_bulk but for edge points on boundary sides.
//    template <unsigned int dim>
//    std::shared_ptr<BoundaryIntegral> add_boundary(const Quadrature &);

    /// Return maximal size of evaluation points objects.
    inline unsigned int max_size() const {
        return max_size_;
    }

    inline void clear() {
        for (uint i=0; i<4; ++i)
            dim_eval_points_[i].clear();
    }

    void create_integrals(std::vector<DimIntegrals> integrals_vec);

private:
    /// Subobject holds evaluation points data of one dimension (0,1,2,3)
    class DimEvalPoints {
    public:
        /// Constructor
        DimEvalPoints(unsigned int dim);

        /// Return size of evaluation points object (number of points).
        inline unsigned int size() const {
            if (dim_==0) return n_subsets_;
        	else return local_points_.size();
        }

        /// Return local coordinates of given local point.
        template<unsigned int dim>
        inline arma::vec::fixed<dim> local_point(unsigned int local_point_idx) const {
            ASSERT_LT(local_point_idx, this->size());
            return local_points_.vec<dim>(local_point_idx);
        }

        /// Return begin index of appropriate subset data.
        inline int subset_begin(unsigned int idx) const {
            ASSERT_LT(idx, n_subsets());
        	return subset_starts_[idx];
        }

        /// Return end index of appropriate subset data.
        inline int subset_end(unsigned int idx) const {
            ASSERT_LT(idx, n_subsets());
        	return subset_starts_[idx+1];
        }

        /// Return number of local points corresponding to subset.
        inline int subset_size(unsigned int idx) const {
            ASSERT_LT(idx, n_subsets());
        	return subset_starts_[idx+1] - subset_starts_[idx];
        }

        /// Return number of subsets.
        inline unsigned int n_subsets() const {
            return n_subsets_;
        }

        /// Adds set of local point to local_points_ (bulk or side).
    	template <unsigned int dim>
        void add_local_points(const Armor::Array<double> & quad_points);

        /// Adds new subset and its end size to subset_starts_ array.
        uint add_subset();

        inline void clear() {
            local_points_.resize(0);
            n_subsets_ = 0;
        }


    private:
        Armor::Array<double> local_points_;                           ///< Local coords of points vector
        std::array<int, EvalPoints::max_subsets+1> subset_starts_;    ///< Indices of subsets data in local_points_ vector, used size is n_subsets_ + 1
        unsigned int n_subsets_;                                      ///< Number of subset
        unsigned int dim_;                                            ///< Dimension of local points
    };

    inline void set_max_size() {
        max_size_ = std::max( std::max( size(0), size(1) ), std::max( size(2), size(3) ) );
    }


    /// Sub objects of dimensions 0,1,2,3
    std::array<DimEvalPoints, 4> dim_eval_points_;

    /// Maps of all BulkIntegrals of dimensions 0,1,2,3
    IntegralHashMap<BulkIntegral> bulk_integrals_;

    /// Maps of all EdgeIntegrals of dimensions 1,2,3
    IntegralHashMap<EdgeIntegral> edge_integrals_;

    /// Maximal number of used EvalPoints.
    unsigned int max_size_;

};


#endif /* EVAL_POINTS_HH_ */
