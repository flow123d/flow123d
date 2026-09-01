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
#include "fem/integral_data.hh"
#include "mesh/range_wrapper.hh"
#include "system/asserts.hh"
#include "system/armor.hh"

class Side;
class Quadrature;
class ElementCacheMap;
class BulkIntegral;
class EdgeIntegral;
class CouplingIntegral;
class BoundaryIntegral;
template <unsigned int dim> class BulkIntegralAcc;
template <unsigned int dim> class EdgeIntegralAcc;
template <unsigned int dim> class CouplingIntegralAcc;
template <unsigned int dim> class BoundaryIntegralAcc;
template <int spacedim> class ElementAccessor;
template <unsigned int spacedim> class PatchFEValues;
namespace internal_integrals {
    class Bulk;
    class Edge;
}


/// Distinguishes quadrature points of DimEvalPoints by type
enum points_domain
{
	bulk_points =0,
	side_points =1,
	repeated_side_points =2
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
     *
     * Returns an object referencing to the EvalPoints and list of its points.
     */
    template <unsigned int dim>
    std::shared_ptr<internal_integrals::Bulk> add_bulk_internal(Quadrature *);

    /// The same as add_bulk but for edge points on sides.
    template <unsigned int dim>
    std::shared_ptr<internal_integrals::Edge> add_edge_internal(Quadrature *);

    /// Return maximal size of evaluation points objects.
    inline unsigned int max_size() const {
        return max_size_;
    }

    inline void clear() {
        for (uint i=0; i<4; ++i)
            dim_eval_points_[i].clear();
    }

    /// Return maximal size of quadrature of given dimension of bulk integral (Bulk, Coupling /lower-dim/)
    uint get_max_bulk_quad_size(unsigned int dim) const;

    /// Return maximal size of quadrature of given dimension of side integral (Edge, Coupling /higher-dim/, Boundary)
    uint get_max_side_quad_size(unsigned int dim) const;

    /// Return vector of bulk quadratures of given dimension from bulk integrals (Bulk, Coupling /lower-dim/)
    std::vector<Quadrature *> get_bulk_quad_vector(unsigned int dim) const;

    /// Return vector of side quadratures of given dimension from side integrals (Edge, Coupling /higher-dim/, Boundary)
    std::vector<Quadrature *> get_side_quad_vector(unsigned int dim) const;

    /// Return begin index of appropriate subset data.
    inline points_domain point_domain(unsigned int dim, unsigned int local_point_idx) const {
    	return dim_eval_points_[dim].point_domain(local_point_idx);
    }

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
        uint add_subset(points_domain point_domain, unsigned int quad_size, unsigned int repeated_points=0);

        inline void clear() {
            local_points_.resize(0);
            n_subsets_ = 0;
        }

        /// Return begin index of appropriate subset data.
        inline points_domain point_domain(unsigned int local_point_idx) const {
            ASSERT_LT(local_point_idx, points_domains_.size());
        	return points_domains_[local_point_idx];
        }

    private:
        Armor::Array<double> local_points_;                           ///< Local coords of points vector
        std::array<int, EvalPoints::max_subsets+1> subset_starts_;    ///< Indices of subsets data in local_points_ vector, used size is n_subsets_ + 1
        std::vector<points_domain> points_domains_;                   ///< Flags hold if quadrature points are bulk or side, temporary data member of FieldFE patch operations
        unsigned int n_subsets_;                                      ///< Number of subset
        unsigned int dim_;                                            ///< Dimension of local points

        friend class EvalPoints;
    };

    inline void set_max_size() {
        max_size_ = std::max( std::max( size(0), size(1) ), std::max( size(2), size(3) ) );
    }

    /// Common implementation of get_max_bulk_quad_size and get_max_side_quad_size
    template<class Integral>
    uint get_max_integral_quad_size(IntegralPtrMap<Integral> integrals, unsigned int dim) const;

    /**
     * Common implementation of get_bulk_quad_vector and get_side_quad_vector
     *
     * Temporary method that allows implementation PatchFeValues operations in FieldFE
     */
    template<class Integral>
    std::vector<Quadrature *> get_quad_vector(IntegralPtrMap<Integral> integrals, unsigned int dim) const;


    /// Sub objects of dimensions 0,1,2,3
    std::array<DimEvalPoints, 4> dim_eval_points_;

    /// Maps of all BulkIntegrals of dimensions 0,1,2,3
    IntegralPtrMap<internal_integrals::Bulk> bulk_integrals_;

    /// Maps of all EdgeIntegrals of dimensions 1,2,3
    IntegralPtrMap<internal_integrals::Edge> edge_integrals_;

    /// Maximal number of used EvalPoints.
    unsigned int max_size_;

};


#endif /* EVAL_POINTS_HH_ */
