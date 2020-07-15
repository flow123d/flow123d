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
#include "mesh/accessors.hh"

class EvalPoints;
class ElementCacheMap;
class DHCellAccessor;
class DHCellSide;


/**
 * @brief Class holds precomputed field values of selected element set.
 *
 * Every field in equation use own instance for every dimension of elements
 * (typically 3 instances for dim = 1,2,3).
 */
template<class elm_type>
class FieldValueCache {
public:
    /// Constructor
    FieldValueCache(unsigned int n_rows, unsigned int n_cols);

    /// Destructor
    ~FieldValueCache();

    /// Reinitialize size of cache
    void reinit(const ElementCacheMap &cache_map);

    /// Return size of data cache (number of stored field values)
    inline unsigned int size() const {
        return data_.size();

    }

    /// Return data vector.
    inline const Armor::Array<elm_type> &data() const
    {
        return data_;
    }

    inline Armor::Array<elm_type> &data()
    {
        return data_;
    }

    /// Return data vector.
    template<uint nr, uint nc = 1>
    typename arma::Mat<elm_type>::template fixed<nr, nc> &get(uint i) {
        return data_.template mat<nr, nc>(i);
    }

    /// Return maximal number of data stored in cache.
    inline unsigned int max_size() const {
        return data_.size();
    }

    /// Return value of evaluation point given by DHCell and local point idx in EvalPoints.
    template<class Value>
    typename Value::return_type get_value(const ElementCacheMap &map,
            const DHCellAccessor &dh_cell, unsigned int eval_points_idx);

    /**
     * Return value of evaluation point given by ElementAccessor and local point idx in EvalPoints.
     *
     * Temporary overload of previous method used on boundary elements.
     */
    template<class Value>
    typename Value::return_type get_value(const ElementCacheMap &map,
            const ElementAccessor<3> elm, unsigned int eval_points_idx);

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

    /// Maximal number of points stored in cache.
    unsigned int n_cache_points_;
};


/**
 * @brief Directing class of FieldValueCache.
 *
 * Manage storing and updating element data (elements of same dimension) to cache. We need only
 * one shared instance of this class for all fields in equation (but typically for dim = 1,2,3).
 *
 * IMPORTANT: Because there are combuned bulk and boundary elements, we must use mesh_idx value
 * to correct identification of elements.
 *
 * TODO: The logic of creating and updating this class is quite complex, describe in which order
 * the methods are supposed to be called and which internal structures are updated when.
 */
class ElementCacheMap {
public:
    /// Number of cached elements which values are stored in cache.
    static constexpr unsigned int n_cached_elements = 20;

    /// Index of invalid element in cache.
    static const unsigned int undef_elem_idx;

    /// Size of block (evaluation of FieldFormula) must be multiple of this value.
    /// TODO We should take this value from BParser and it should be dependent on processor configuration.
    static const unsigned int simd_size_double;

    /**
     * Holds elements indices of one region stored in cache.
     *
     * TODO: Auxilliary structure, needs optimization:
     * Proposal -
     *  - store element indices of all regions to one unique set (hold as data member of ElementCacheMap)
     *  - sort elements by regions in prepare_elements_to_update method
     *  - there is only problem that we have to construct ElementAccessors repeatedly
     */
    struct RegionData {
    public:
        /// Constructor
        RegionData() : n_elements_(0), cache_position_(ElementCacheMap::undef_elem_idx) {}

        /// Add element if does not exist
        bool add(ElementAccessor<3> elm) {
            if (std::find(elm_indices_.begin(), elm_indices_.begin()+n_elements_, elm.mesh_idx()) == elm_indices_.begin()+n_elements_) {
                elm_indices_[n_elements_] = elm.mesh_idx();
                n_elements_++;
                return true;
            } else
                return false;
        }

        /// Array of elements mesh_idx, ensures element uniqueness
        std::array<unsigned int, ElementCacheMap::n_cached_elements> elm_indices_;
        /// Number of element indices
        unsigned int n_elements_;
        /// Holds positions of regions in cache
		/// Elements of the common region forms continuous chunks in the
        /// ElementCacheMap table. This array gives start indices of the regions
        /// in array of all cached elements.
        unsigned int cache_position_;
    };

    /**
     * Holds data of regions (and their elements) stored in ElementCacheMap.
     *
     * TODO: Needs further optimization.
     * Is this like a public par of the data?
     * The name is too generic.
     */
    class UpdateCacheHelper {
    public:
        /// Maps of data of different regions in cache
    	/// TODO: auxiliary data membershould be removed or moved to sepaate structure.
        std::unordered_map<unsigned int, RegionData> region_cache_indices_map_;

        /// Maps of begin and end positions of different regions data in FieldValueCache
        std::array<unsigned int, ElementCacheMap::n_cached_elements+1> region_value_cache_range_;

        /// Maps of begin and end positions of elements of different regions in ElementCacheMap
        std::array<unsigned int, ElementCacheMap::n_cached_elements+1> region_element_cache_range_;
    };

    /// Constructor
    ElementCacheMap();

    /// Destructor
    ~ElementCacheMap();

    /// Init cache
    void init(std::shared_ptr<EvalPoints> eval_points);

    /// Adds element to region_cache_indices_map_ set.
    void add(const DHCellAccessor &dh_cell);

    /// Same as previous but passes DHCellSide as parameter.
    void add(const DHCellSide &cell_side);

    /// Adds element to region_cache_indices_map_ set.
    void add(const ElementAccessor<3> &elm_acc);

    /// Prepare data member before reading data to cache.
    void prepare_elements_to_update();

    /// Create map of used eval points on cached elements.
    void create_elements_points_map();

    /// Start update of cache.
    void start_elements_update();

    /// Finish update after reading data to cache.
    void finish_elements_update();

    /// Return update cache data helper
    inline const UpdateCacheHelper &update_cache_data() const {
        return update_data_;
    }

    /// Return update cache data helper
    inline UpdateCacheHelper &update_cache_data() {
        return update_data_;
    }

    /// Getter of eval_points object.
    inline std::shared_ptr<EvalPoints> eval_points() const {
        return eval_points_;
    }

    /**
     * Set used element eval points to value ElementCacheMap::point_in_proggress
     *
     * @param dh_cell     Specified element
     * @param subset_idx  Index of subset
     * @param data_size   Number of points
     * @param start_point Index of first used point in subset (e.g. subset holds eval points of all sides but EdgeIntegral represents only one of them)
     */
    void mark_used_eval_points(const DHCellAccessor &dh_cell, unsigned int subset_idx, unsigned int data_size, unsigned int start_point=0);

    /**
     * Same as previous but passes ElementAccessor.
     * Temporary method used on boundary elements.
     */
    void mark_used_eval_points(const ElementAccessor<3> elm, unsigned int subset_idx, unsigned int data_size, unsigned int start_point=0);

    /*
     * Return index of point in FieldValueCache.
     *
     * @param cache_elm_idx  idx of ElementAccessor in ElementCacheMap
     * @param loc_point_idx  Index of local point in EvalPoints
     */
    inline int get_field_value_cache_index(unsigned int cache_elm_idx, unsigned int loc_point_idx) const {
        ASSERT_PTR_DBG(element_eval_points_map_);
        return element_eval_points_map_[cache_elm_idx][loc_point_idx];
    }

    /// Return mesh_idx of element stored at given position of ElementCacheMap
    inline unsigned int elm_idx_on_position(unsigned pos) const {
        return elm_idx_[pos];
    }

    /// Return position of element stored in ElementCacheMap
    inline unsigned int position_in_cache(unsigned mesh_elm_idx) const {
        std::unordered_map<unsigned int, unsigned int>::const_iterator it = cache_idx_.find(mesh_elm_idx);
        if ( it != cache_idx_.end() ) return it->second;
        else return ElementCacheMap::undef_elem_idx;
    }

    /// Return number of stored elements.
    inline unsigned int n_elements() const {
        return update_data_.region_element_cache_range_[update_data_.region_cache_indices_map_.size()];
    }

    inline unsigned int region_chunk(unsigned int region_idx) const {
        return update_data_.region_cache_indices_map_.find(region_idx)->second.cache_position_;
    }

    /// Set index of cell in ElementCacheMap (or undef value if cell is not stored in cache).
    DHCellAccessor & operator() (DHCellAccessor &dh_cell) const;
protected:

    /// Special constant (@see element_eval_points_map_).
    static const int unused_point = -2;

    /// Special constant (@see element_eval_points_map_).
    static const int point_in_proggress = -1;

    /// Reset all items of elements_eval_points_map
    void clear_element_eval_points_map();

    /// Add element to appropriate region data of update_data_ object
    void add_to_region(ElementAccessor<3> elm);

    /// Vector of element indexes stored in cache.
    /// TODO: could be moved to UpdateCacheHelper structure
    std::vector<unsigned int> elm_idx_;

    /// Map of element indices stored in cache, allows reverse search to previous vector.
    /// TODO: could be moved to UpdateCacheHelper structure
    std::unordered_map<unsigned int, unsigned int> cache_idx_;

    /// Pointer to EvalPoints
    std::shared_ptr<EvalPoints> eval_points_;

    /// Holds data used for cache update.
    UpdateCacheHelper update_data_;

    /// Flag is set down during update of cache when this can't be read
    bool ready_to_reading_;

    /**
     * Two dimensional array provides indexes to FieldValueCache.
     *
     * 1: Over elements holds in ElementCacheMap
     * 2: Over EvalPoints for each element
     *
     * Array is filled in those three steps:
     * a. Reset - all items are set to ElementCacheMap::unused_point
     * b. Used eval points are set to ElementCacheMap::point_in_proggress
     * c. Eval points marked in previous step are sequentially numbered
     */
    // TODO: What are the dimensions of the table?
    // should be n_cached_elements * n_eval_points, document it.
    //
    // Better use just int *, and use just single allocation of the whole table
    // current impl. have bad memory locality. Define a private access method.
    int **element_eval_points_map_;
};



#endif /* FIELD_VALUE_CACHE_HH_ */
