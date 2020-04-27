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
#include "mesh/range_wrapper.hh"
#include "system/asserts.hh"
#include "system/armor.hh"

class Side;
class Quadrature;
class BulkIntegral;
class EdgeIntegral;
class CouplingIntegral;
class BoundaryIntegral;
template <int spacedim> class ElementAccessor;


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
        return dim_eval_points_[dim-1].size();
    }

    /// Return local coordinates of given local point and appropriate dim.
    template<unsigned int dim>
    inline arma::vec::fixed<dim> local_point(unsigned int local_point_idx) const {
        return dim_eval_points_[dim-1].local_point<dim>(local_point_idx);
    }

    /// Return begin index of appropriate subset data.
    inline int subset_begin(unsigned int dim, unsigned int idx) const {
        return dim_eval_points_[dim-1].subset_begin(idx);
    }

    /// Return end index of appropriate subset data.
    inline int subset_end(unsigned int dim, unsigned int idx) const {
        return dim_eval_points_[dim-1].subset_end(idx);
    }

    /// Return number of local points corresponding to subset.
    inline int subset_size(unsigned int dim, unsigned int idx) const {
        return dim_eval_points_[dim-1].subset_size(idx);
    }

    /// Return number of subsets.
    inline unsigned int n_subsets(unsigned int dim) const {
        return dim_eval_points_[dim-1].n_subsets();
    }

    /**
     * Registers point set from quadrature.
     * Returns an object referencing to the EvalPoints and list of its points.
     */
    template <unsigned int dim>
    std::shared_ptr<BulkIntegral> add_bulk(const Quadrature &);

    /// The same as add_bulk but for edge points on sides.
    template <unsigned int dim>
    std::shared_ptr<EdgeIntegral> add_edge(const Quadrature &);

    /// The same as add_bulk but for points between side points of element of dim and bulk points of element of dim-1.
    template <unsigned int dim>
    std::shared_ptr<CouplingIntegral> add_coupling(const Quadrature &);

    /// The same as add_bulk but for edge points on boundary sides.
    template <unsigned int dim>
    std::shared_ptr<BoundaryIntegral> add_boundary(const Quadrature &);

    /// Return maximal size of evaluation points objects .
    inline unsigned int max_size() const {
        return std::max( size(1), std::max( size(2), size(3) ) );
    }

private:
    class DimEvalPoints {
    public:
        /// Constructor
        DimEvalPoints(unsigned int dim);

        /// Return size of evaluation points object (number of points).
        inline unsigned int size() const {
            return local_points_.size();
        }

        /// Return local coordinates of given local point.
        template<unsigned int dim>
        inline arma::vec::fixed<dim> local_point(unsigned int local_point_idx) const {
            ASSERT_LT_DBG(local_point_idx, this->size());
            return local_points_.vec<dim>(local_point_idx);
        }

        /// Return begin index of appropriate subset data.
        inline int subset_begin(unsigned int idx) const {
            ASSERT_LT_DBG(idx, n_subsets());
        	return subset_starts_[idx];
        }

        /// Return end index of appropriate subset data.
        inline int subset_end(unsigned int idx) const {
            ASSERT_LT_DBG(idx, n_subsets());
        	return subset_starts_[idx+1];
        }

        /// Return number of local points corresponding to subset.
        inline int subset_size(unsigned int idx) const {
            ASSERT_LT_DBG(idx, n_subsets());
        	return subset_starts_[idx+1] - subset_starts_[idx];
        }

        /// Return number of subsets.
        inline unsigned int n_subsets() const {
            return n_subsets_;
        }

        /// Adds set of local point to local_points_ (bulk or side of given permutation).
    	template <unsigned int dim>
        void add_local_points(const Armor::Array<double> & quad_points);

        /// Find position of local point (coords) in subvector of local points given by limits <data_begin,  ... data_end)
    	template <unsigned int dim>
        unsigned int find_permute_point(arma::vec coords, unsigned int data_begin, unsigned int data_end);

        /// Adds new subset and its end size to subset_starts_ array.
        void add_subset();
    private:
        Armor::Array<double> local_points_;                           ///< Local coords of points vector
        std::array<int, EvalPoints::max_subsets+1> subset_starts_;    ///< Indices of subsets data in local_points_ vector, used size is n_subsets_ + 1
        unsigned int n_subsets_;                                      ///< Number of subset
        unsigned int dim_;                                            ///< Dimension of local points
    };

    /// Sub objects of dimensions 1,2,3
    std::array<DimEvalPoints, 3> dim_eval_points_;

    friend class EvalSubSet;
};


#endif /* EVAL_POINTS_HH_ */
