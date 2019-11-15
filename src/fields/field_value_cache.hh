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
class DHCellAccessor;


template<class Value>
class FieldValueCache {
public:
	typedef typename Value::element_type elm_type;

    /// Default constructor
    FieldValueCache();

    /// Constructor
    FieldValueCache(const EvalPoints *eval_points, unsigned int n_cache_points);

    /// Destructor
    ~FieldValueCache();

    /// Marks the used local points
    void mark_used(EvalSubset sub_set);

    /// Getter for used_points
    inline const std::set<int> &used_points() const {
        return used_points_;
    }

    /// Return size of data cache (number of stored field values)
    inline unsigned int size() const {
        return data_.n_vals();

    }

    /// Return data vector.
    inline Armor::Array<elm_type> &data() {
        return data_;
    }

    /// Return data vector.
    inline const EvalPoints *eval_points() const {
        return eval_points_;
    }

private:
    /// Data cache
    Armor::Array<elm_type> data_;

    /// Holds indices of used local points
    std::set<int> used_points_;  // TODO: test unorderd_set during tuning of the performance

    /// Pointer to EvalPoints
    const EvalPoints *eval_points_;

    /// Dimension (control data member)
    unsigned int dim_;
};


class ElementCacheMap {
public:
    /// Number of cached elements which values are stored in cache.
    static const unsigned int n_cached_elements;

    /// Index of invalid element in cache.
    static const unsigned int undef_elem_idx;

    /// Constructor
    ElementCacheMap(unsigned int dim);

    /// Adds element to added_elements_ set.
    void add(DHCellAccessor dh_cell);

    /// Clean helper data member after reading data to cache.
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
