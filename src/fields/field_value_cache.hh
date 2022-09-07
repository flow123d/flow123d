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
#include "tools/mixed.hh"
#include "tools/revertable_list.hh"
#include "fem/dofhandler.hh"

class EvalPoints;
class ElementCacheMap;
class DHCellAccessor;
class DHCellSide;
template < template<IntDim...> class DimAssembly> class GenericAssembly;
template < template<IntDim...> class DimAssembly> class GenericAssemblyObserve;


/**
 * @brief Class holds precomputed field values of selected element set.
 *
 * Every field in equation use own instance used for elements of all dimensions.
 */
template<class elm_type> using FieldValueCache = Armor::Array<elm_type>;


/**
 * Specifies eval points by idx of region, element and eval point.
 *
 * TODO Add better description after finish implementation
 */
struct EvalPointData {
    EvalPointData() {}              ///< Default constructor
    /// Constructor sets all data members
    EvalPointData(unsigned int i_reg, unsigned int i_ele, unsigned int i_ep, unsigned int dh_loc_idx)
    : i_reg_(i_reg), i_element_(i_ele), i_eval_point_(i_ep), dh_loc_idx_(dh_loc_idx) {}
    /// Copy constructor
    EvalPointData(const EvalPointData &other)
    : i_reg_(other.i_reg_), i_element_(other.i_element_), i_eval_point_(other.i_eval_point_), dh_loc_idx_(other.dh_loc_idx_) {}


    bool operator < (const EvalPointData &other) {
        if (i_reg_ == other.i_reg_) {
            if (i_element_ == other.i_element_)
                return (i_eval_point_ < other.i_eval_point_);
            else
                return (i_element_ < other.i_element_);
        } else
            return (i_reg_ < other.i_reg_);
    }

    unsigned int i_reg_;            ///< region_idx of element
    unsigned int i_element_;        ///< mesh_idx of ElementAccessor appropriate to element
    unsigned int i_eval_point_;     ///< index of point in EvalPoint object
    unsigned int dh_loc_idx_;       ///< local index of cell in DOF handler
};




/**
 * @brief Auxiliary data class holds number of elements in cache and allow to set this value
 * explicitly (e.g. as input parameter).
 *
 * Implementation is done as singletone with two access through static methods 'get' and 'set'.
 *
 * TODO: This is actually a constant so make it a constant int in element cache map.
 */
class CacheMapElementNumber {
public:
	/// Return number of stored elements
	static unsigned int get() {
	    return get_instance().n_elem_;
	}

	/// Set number of stored elements
	static void set(unsigned int n_elem) {
	    get_instance().n_elem_ = n_elem;
	}

	CacheMapElementNumber(CacheMapElementNumber const&) = delete;  ///< We don't need copy constructor.
	void operator=(CacheMapElementNumber const&)        = delete;  ///< We don't need assignment operator.

private:
	/// Forbiden default constructor
	CacheMapElementNumber() : n_elem_(300) {}


    static CacheMapElementNumber& get_instance()
    {
        static CacheMapElementNumber instance;
        return instance;
    }

    /// Maximal number of elements stored in cache.
    unsigned int n_elem_;
};


/**
 * @brief Directing class of FieldValueCache.
 *
 * Manage storing and updating element data (elements of same dimension) to cache. We need only
 * one shared instance of this class for all fields in equation (but typically for dim = 1,2,3).
 *
 * IMPORTANT: Because there are combined bulk and boundary elements, we must use mesh_idx value
 * to correct identification of elements.
 *
 * TODO:
 * 1. Generic assembly pass through the patch collecting needed quadrature points. (PASS ORDER)
 * 2. Then we sort these points for efficient chae_upadate of the fields (CACHE ORDER)
 * 3. We pass through the patch again evaluating actual integrals. This second pass is currently inefficient since
 *    we can not map efficiently from the PASS ORDER to the CACHE ORDER), this leads to many complications in
 *    quad point classes.
 * We should:
 * 1. have templated patch iteration mechanism, so that we can iterate through the evaluated integrals
 *    twice in consistent way performing:
 *    FIRST: collection of evaluation points
 *    SECOND: evaluation of integrals using the fields and fe_values
 *    Having a consistent implementation allows us to assign unique indices to the integral points on the patch.
 * 2. Sort collected points, remove duplicities, mark new indices to the original list of points.
 * 3. SECOND pass, use unique index to find the point in the cache.
 *
 * Resulting simplifications:
 * - no need for associating point operations for Edge, Coupling and Boundary points,
 *   we add  assiciated eval point pairs as separate eval points with unique ids. Consistent integral iteration
 *   allows us to simply take two succesive points at the second pass.
 * - no need for the matrix mapping (element, eval_point)  to the cache index
 */
class ElementCacheMap {
public:
    /// Index of invalid element in cache.
    static const unsigned int undef_elem_idx;

    /// Size of block (evaluation of FieldFormula) must be multiple of this value.
    /// TODO We should take this value from BParser and it should be dependent on processor configuration.
    static const unsigned int simd_size_double;

    /// Constructor
    ElementCacheMap();

    /// Destructor
    ~ElementCacheMap();

    /// Init cache
    void init(std::shared_ptr<EvalPoints> eval_points);

    /// Create patch of cached elements before reading data to cache.
    void create_patch();

    /// Reset all items of elements_eval_points_map
    inline void clear_element_eval_points_map() {
        ASSERT_PTR(element_eval_points_map_);
        unsigned int last_element_idx = -1, i_elem_row = -1;
        for (unsigned int i=0; i<eval_point_data_.permanent_size(); ++i) {
            if (eval_point_data_[i].i_element_ != last_element_idx) { // new element
                i_elem_row++;
            	last_element_idx =eval_point_data_[i].i_element_;
            }
        	set_element_eval_point(i_elem_row, eval_point_data_[i].i_eval_point_, ElementCacheMap::unused_point);
        }
        eval_point_data_.reset();
    }

    /// Start update of cache.
    void start_elements_update();

    /// Finish update after reading data to cache.
    void finish_elements_update();

    /// Getter of eval_points object.
    inline std::shared_ptr<EvalPoints> eval_points() const {
        return eval_points_;
    }

    /** Adds EvalPointData using emplace_back.
     *  Arguments correspond to constructor of EvalPointData.
     */
    inline void add_eval_point(unsigned int i_reg, unsigned int i_ele, unsigned int i_eval_point, unsigned int dh_loc_idx)
    {
        eval_point_data_.emplace_back(i_reg, i_ele, i_eval_point, dh_loc_idx);
        set_of_regions_.insert(i_reg);
    }

    /// Returns number of eval. points with addition of max simd duplicates due to regions. 
    inline unsigned int get_simd_rounded_size()
    {
        return eval_point_data_.temporary_size() + (simd_size_double - 1)*set_of_regions_.size();
    }

    /*
     * Access to item of \p element_eval_points_map_ like to two-dimensional array.
     *
     * @param i_elem_in_cache  idx of ElementAccessor in ElementCacheMap
     * @param i_eval_point     index of local point in EvalPoints
     * @return                 index of point in FieldValueCache.
     */
    inline int element_eval_point(unsigned int i_elem_in_cache, unsigned int i_eval_point) const {
        ASSERT_PTR(element_eval_points_map_);
        return element_eval_points_map_[i_elem_in_cache*eval_points_->max_size()+i_eval_point];
    }

    /// Return mesh_idx of element stored at given position of ElementCacheMap
    inline unsigned int elm_idx_on_position(unsigned pos) const {
        return elm_idx_[pos];
    }

    /// Return position of element stored in ElementCacheMap
    inline unsigned int position_in_cache(unsigned mesh_elm_idx, bool bdr=false) const {
        std::unordered_map<unsigned int, unsigned int>::const_iterator it;
        if (bdr) {
            it = element_to_map_bdr_.find(mesh_elm_idx);
            if ( it != element_to_map_bdr_.end() ) return it->second;
            else return ElementCacheMap::undef_elem_idx;
        } else {
            it = element_to_map_.find(mesh_elm_idx);
            if ( it != element_to_map_.end() ) return it->second;
            else return ElementCacheMap::undef_elem_idx;
        }
    }

    /// Return number of stored regions.
    inline unsigned int n_regions() const {
        return regions_starts_.permanent_size() - 1;
    }

    /// Return number of stored elements.
    inline unsigned int n_elements() const {
        return element_starts_.permanent_size() - 1;
    }

    /// Return begin position of element chunk in FieldValueCache
    inline unsigned int element_chunk_begin(unsigned int elm_patch_idx) const {
        ASSERT_LT(elm_patch_idx, n_elements());
        return element_starts_[elm_patch_idx];
    }

    /// Return end position of element chunk in FieldValueCache
    inline unsigned int element_chunk_end(unsigned int elm_patch_idx) const {
        ASSERT_LT(elm_patch_idx, n_elements());
        return element_starts_[elm_patch_idx+1];
    }

    /// Return begin position of region chunk in FieldValueCache
    inline unsigned int region_chunk_begin(unsigned int region_patch_idx) const {
        ASSERT_LT(region_patch_idx, n_regions());
        return element_starts_[ regions_starts_[region_patch_idx] ];
    }

    /// Return end position of region chunk in FieldValueCache
    inline unsigned int region_chunk_end(unsigned int region_patch_idx) const {
        ASSERT_LT(region_patch_idx, n_regions());
        return element_starts_[ regions_starts_[region_patch_idx+1] ];
    }

    /// Return begin position of region chunk specified by position in map
    inline unsigned int region_chunk_by_map_index(unsigned int r_idx) const {
        if (r_idx <= n_regions()) return element_starts_[ regions_starts_[r_idx] ];
        else return ElementCacheMap::undef_elem_idx;
    }

    /// Return begin position of region chunk specified by position in map
    inline unsigned int region_idx_from_chunk_position(unsigned int chunk_pos) const {
    	return eval_point_data_[ this->region_chunk_by_map_index(chunk_pos) ].i_reg_;
    }

    /// Return item of eval_point_data_ specified by its position
    inline const EvalPointData &eval_point_data(unsigned int point_idx) const {
        return eval_point_data_[point_idx];
    }

    /// Return value of evaluation point given by idx of element in patch and local point idx in EvalPoints from cache.
    template<class Value>
    inline typename Value::return_type get_value(const FieldValueCache<typename Value::element_type> &field_cache,
            unsigned int elem_patch_idx, unsigned int eval_points_idx) const {
        ASSERT_EQ(Value::NRows_, field_cache.n_rows());
        ASSERT_EQ(Value::NCols_, field_cache.n_cols());
        unsigned int value_cache_idx = this->element_eval_point(elem_patch_idx, eval_points_idx);
        ASSERT(value_cache_idx != ElementCacheMap::undef_elem_idx);
        return Value::get_from_array(field_cache, value_cache_idx);
    }
protected:

    /// Special constant (@see element_eval_points_map_).
    static const int unused_point = -1;

    /// Base number of stored regions in patch
    static const unsigned int regions_in_chunk = 3;

    /// Base number of stored elements in patch
    static const unsigned int elements_in_chunk = 10;

    /// Set item of \p element_eval_points_map_.
    inline void set_element_eval_point(unsigned int i_elem_in_cache, unsigned int i_eval_point, int val) const {
        ASSERT_PTR(element_eval_points_map_);
        element_eval_points_map_[i_elem_in_cache*eval_points_->max_size()+i_eval_point] = val;
    }

    /// Vector of element indexes stored in cache.
    std::vector<unsigned int> elm_idx_;

    /// Pointer to EvalPoints
    std::shared_ptr<EvalPoints> eval_points_;

    /// Flag is set down during update of cache when this can't be read
    bool ready_to_reading_;

    /**
     * This array provides indexes to FieldValueCache.
     *
     * This one dimensional array behaves like two dimensional factually.
     * Size is set to 'n_cached_elements * n_eval_points' and items are
     * accessible through two indices:
     *
     * 1: Over elements holds in ElementCacheMap
     * 2: Over EvalPoints for each element
     *
     * Use always and only methods \p element_eval_point for read and
     * \p set_element_eval_point (for write) to access to items!
     *
     * Array is filled in those three steps:
     * a. Reset - all items are set to ElementCacheMap::unused_point
     * b. Used eval points are set to ElementCacheMap::point_in_proggress
     * c. Eval points marked in previous step are sequentially numbered
     *
     * TODO improve description
     */
    int *element_eval_points_map_;

    ///< Holds data of evaluating points in patch.
    RevertableList<EvalPointData> eval_point_data_;

    /// @name Holds start positions and orders of region chunks and element chunks
    // @{

    RevertableList<unsigned int> regions_starts_;         ///< Start positions of elements in regions (size = n_regions+1, last value is end of last region)
    RevertableList<unsigned int> element_starts_;         ///< Start positions of elements in eval_point_data_ (size = n_elements+1)
    std::unordered_map<unsigned int, unsigned int> element_to_map_;     ///< Maps bulk element_idx to element index in patch - TODO remove
    std::unordered_map<unsigned int, unsigned int> element_to_map_bdr_; ///< Maps boundary element_idx to element index in patch - TODO remove

    // @}

    /// Keeps set of unique region indices of added eval. points.
    std::unordered_set<unsigned int> set_of_regions_;

    // TODO: remove friend class
    template < template<IntDim...> class DimAssembly>
    friend class GenericAssembly;
    template < template<IntDim...> class DimAssembly>
    friend class GenericAssemblyObserve;
};



#endif /* FIELD_VALUE_CACHE_HH_ */
