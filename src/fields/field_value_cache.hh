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
 * @file    field_value_cache.hh
 * @brief
 * @author  David Flanderka
 */

#ifndef FIELD_VALUE_CACHE_HH_
#define FIELD_VALUE_CACHE_HH_

#include <set>
#include <unordered_map>
#include <unordered_set>
#include <vector>
#include "system/armor.hh"
#include "fields/eval_points.hh"

class EvalSubset;
class DHCellAccessor;
class Mesh;
template <int spacedim> class ElementAccessor;


/**
 * @brief Class holds precomputed field values of selected element set.
 *
 * Every field in equation use own instance for every dimension of elements
 * (typically 3 instances for dim = 1,2,3).
 */
template<class elm_type, class return_type>
class FieldValueCache {
public:
    /// Constructor
    FieldValueCache(unsigned int n_rows, unsigned int n_cols);

    /// Destructor
    ~FieldValueCache();

    /// Initialize cache
    void init(std::shared_ptr<EvalPoints> eval_points, unsigned int dim, unsigned int n_cache_elements);

    /// Marks the used local points
    void mark_used(std::shared_ptr<EvalSubset> sub_set);

    /// Getter for subsets of used points
    inline const std::array<bool, EvalPoints::max_subsets> &used_subsets() const {
        return used_subsets_;
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
        return eval_points_->n_subsets(dim_);
    }

    /// Return size of data cache (number of stored field values)
    inline unsigned int size() const {
        return data_.size();

    }

    /// Return dimension
    inline unsigned int dim() const {
        return dim_;
    }

    /// Return data vector.
    inline Armor::Array<elm_type> &data() {
        return data_;
    }

    /// Return data vector.
    template<uint nr, uint nc = 1>
    inline typename arma::Mat<elm_type>::template fixed<nr, nc> &get(uint i) {
        return data_.mat<nr, nc>(i);
    }

    /// Return data vector.
    inline std::shared_ptr<EvalPoints> eval_points() const {
        return eval_points_;
    }

    /// Return number of elements that data is stored in cache.
    inline unsigned int n_cache_elements() const {
        return n_cache_elements_;
    }

    /// Return value of evaluation point given by DHCell and local point idx in EvalPoints.
    template<uint nRows, uint nCols>
    typename arma::Mat<elm_type>::template fixed<nRows, nCols> get_value(DHCellAccessor dh_cell, unsigned int subset_idx, unsigned int eval_points_idx);

private:
    /**
     * Data cache.
     *
     * Data is ordered like three dimensional table. The highest level is determinated by subsets,
     * those data ranges are holds in subset_starts. Data block size of each subset is determined
     * by number of eval_points (of subset) and maximal number of stored elements.
     * The table is allocated to hold all subsets, but only those marked in used_subsets are updated.
     * Order of subsets is same as in eval_points.
     */
    Armor::Array<elm_type> data_;

    /// Holds if blocks of local points are used or not.
    std::array<bool, EvalPoints::max_subsets> used_subsets_;

    /// Indices of subsets begin and end position in data array.
    std::array<int, EvalPoints::max_subsets+1> subset_starts_;

    /// Pointer to EvalPoints.
    std::shared_ptr<EvalPoints> eval_points_;

    /// Maximal number of elements stored in cache.
    unsigned int n_cache_elements_;

    /// Dimension (control data member).
    unsigned int dim_;
};


/**
 * @brief Directing class of FieldValueCache.
 *
 * Manage storing and updating element data (elements of same dimension) to cache. We need only
 * one shared instance of this class for all fields in equation (but typically for dim = 1,2,3).
 */
class ElementCacheMap {
public:
	typedef std::vector<ElementAccessor<3>> ElementSet;

    /// Number of cached elements which values are stored in cache.
    static const unsigned int n_cached_elements;

    /// Index of invalid element in cache.
    static const unsigned int undef_elem_idx;

    /// Holds helper data necessary for cache update.
    class UpdateCacheHelper {
    public:
        /// Set of element indexes that wait for storing to cache.
        std::unordered_set<unsigned int> added_elements_;

        /// Holds old and new indexes of elements preserved in cache.
        std::unordered_map<unsigned int, unsigned int> preserved_elements_;

        /// Maps elements of different regions
        std::unordered_map<unsigned int, ElementSet> region_element_map_;

        /// Maps of begin positions in cache of different regions
        std::unordered_map<unsigned int, unsigned int> region_cache_begin_;
    };

    /// Constructor
    ElementCacheMap();

    /// Init dimension data member
    inline void init(unsigned int dim) {
    	ASSERT_EQ(dim_, EvalPoints::undefined_dim).error("Repeated initialization!");
    	this->ready_to_reading_ = true;
    	this->dim_ = dim;
    }

    /// Adds element to added_elements_ set.
    void add(const DHCellAccessor &dh_cell);

    /// Prepare data member before reading data to cache.
    void prepare_elements_to_update(Mesh *mesh);

    /// Clean helper data after reading data to cache.
    void clear_elements_to_update();

    /// Return dimension
    inline unsigned int dim() const {
        return dim_;
    }

    /// Return update cache data helper
    inline const UpdateCacheHelper &update_cache_data() const {
        return update_data_;
    }

    /// Set index of cell in ElementCacheMap (or undef value if cell is not stored in cache).
    DHCellAccessor & operator() (DHCellAccessor &dh_cell) const;
private:
    /// Vector of element indexes stored in cache.
    std::vector<unsigned int> elm_idx_;

    /// Map of element indices stored in cache, allows reverse search to previous vector.
    std::unordered_map<unsigned int, unsigned int> cache_idx_;

    /// Dimension (control data member)
    unsigned int dim_;

    /// Holds data used for cache update.
    UpdateCacheHelper update_data_;

    /// Flag is set down during update of cache when this can't be read
    bool ready_to_reading_;
};



#endif /* FIELD_VALUE_CACHE_HH_ */
