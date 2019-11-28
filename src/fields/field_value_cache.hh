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

class EvalPoints;
class EvalSubset;
class ElementCacheMap;
class DHCellAccessor;


/**
 * @brief Class holds precomputed field values of selected element set.
 *
 * Every field in equation use own instance for every dimension of elements
 * (typically 3 instances for dim = 1,2,3).
 */
template<class elm_type, class Value>
class FieldValueCache {
private:
    /// Maximal nuber of subsets.
    static const unsigned int max_subsets = 10;
public:
    /// Constructor
    FieldValueCache(unsigned int n_rows, unsigned int n_cols);

    /// Destructor
    ~FieldValueCache();

    /// Constructor
    void init(std::shared_ptr<EvalPoints> eval_points, const ElementCacheMap *cache_map, unsigned int n_cache_points);

    /// Marks the used local points
    void mark_used(std::shared_ptr<EvalSubset> sub_set);

    /// Getter for subsets of used points
    inline const std::array<int, FieldValueCache::max_subsets> &used_subsets() const {
        return used_subsets_;
    }

    /// Return size of data cache (number of stored field values)
    inline unsigned int size() const {
        return data_.n_vals();

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
    inline Armor::Mat<elm_type, nr, nc> &get(uint i) {
        return data_.get<nr, nc>(i);
    }

    /// Return data vector.
    inline std::shared_ptr<EvalPoints> eval_points() const {
        return eval_points_;
    }

private:
    /// Data cache.
    Armor::Array<elm_type> data_;

    /// Holds indices of used blocks of local points.
    std::array<int, FieldValueCache::max_subsets> used_subsets_;

    /// Indices of subsets begin and end position in data array.
    std::array<int, FieldValueCache::max_subsets+1> subset_starts_;

    /// Pointer to EvalPoints.
    std::shared_ptr<EvalPoints> eval_points_;

    /// Pointer to ElementCacheMap.
    const ElementCacheMap *element_cache_map_;

    /// Maximal number of elements stored in cache.
    unsigned int n_cache_points_;

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
    /// Number of cached elements which values are stored in cache.
    static const unsigned int n_cached_elements;

    /// Index of invalid element in cache.
    static const unsigned int undef_elem_idx;

    /// Constructor
    ElementCacheMap();

    /// Init dimension data member
    inline void init(unsigned int dim) {
        this->dim_ = dim;
    }

    /// Adds element to added_elements_ set.
    void add(const DHCellAccessor &dh_cell);

    /// Clean helper data member before reading data to cache.
    void prepare_elements_to_update();

    /// Clean helper data member after reading data to cache.
    void clear_elements_to_update();

    /// Getter for begin_idx_
    inline unsigned int begin_idx() const {
        return begin_idx_;
    }

    /// Getter for end_idx_
    inline unsigned int end_idx() const {
        return end_idx_;
    }

    /// Getter for added_elements_
    inline const std::unordered_set<unsigned int> &added_elements() const {
        return added_elements_;
    }

    /// Return dimension
    inline unsigned int dim() const {
        return dim_;
    }

    /// Set index of cell in ElementCacheMap (or undef value if cell is not stored in cache).
    DHCellAccessor & operator() (DHCellAccessor &dh_cell) const;
private:
    /// Vector of element indexes stored in cache.
    std::vector<unsigned int> elm_idx_;

    /// Map of element indices stored in cache, allows reverse search to previous vector.
    std::unordered_map<unsigned int, unsigned int> cache_idx_;

    /// Set of element indexes that wait for storing to cache.
    std::unordered_set<unsigned int> added_elements_;

    /// Holds index to elm_idx_ vector corresponding to begin index stored in added_elements_ vector.
    unsigned int begin_idx_;

    /// Holds index to elm_idx_ vector corresponding to end index stored in added_elements_ vector.
    unsigned int end_idx_;

    /// Dimension (control data member)
    unsigned int dim_;
};



#endif /* FIELD_VALUE_CACHE_HH_ */
