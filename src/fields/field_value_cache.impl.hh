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
 * @file    field_value_cache.impl.hh
 * @brief
 * @author  David Flanderka
 */

#ifndef FIELD_VALUE_CACHE_IMPL_HH_
#define FIELD_VALUE_CACHE_IMPL_HH_

#include "fields/field_value_cache.hh"
#include "fem/dh_cell_accessor.hh"
#include "mesh/accessors.hh"
#include "mesh/ref_element.hh"
#include "fem/mapping_p1.hh"

template<class elm_type, class return_type>
template<uint nRows, uint nCols>
Armor::Mat<elm_type, nRows, nCols> FieldValueCache<elm_type, return_type>::get_value(DHCellAccessor dh_cell,
        unsigned int subset_idx, unsigned int eval_points_idx) {

    ASSERT(dh_cell.element_cache_index() != ElementCacheMap::undef_elem_idx)(dh_cell.elm_idx());
    unsigned int points_per_element = this->subset_size(subset_idx) / n_cache_elements_;
    unsigned int subset_point_idx = eval_points_idx - eval_points_->subset_begin(subset_idx);
    return data_.get<nRows, nCols>(this->subset_begin(subset_idx) + dh_cell.element_cache_index() * points_per_element + subset_point_idx);
}


template<unsigned int elemdim>
void ElementCacheMap::compute_global_coords(std::shared_ptr<EvalPoints> eval_points, MappingP1<elemdim,3> &mapping, Mesh *mesh) {
    ASSERT_EQ(elemdim, this->dim_);
    if (global_coords_.n_vals() == 0) { // initialization in first using
        global_coords_.resize(eval_points->size()*ElementCacheMap::n_cached_elements);
    } else { //check size
        ASSERT_EQ(global_coords_.n_vals(), eval_points->size()*ElementCacheMap::n_cached_elements);
    }

    for (unsigned int i_elm=0; i_elm<cache_idx_.size(); ++i_elm) { // number of 'active' elements in elm_idx_ is equal to size of cache_idx_
        ElementAccessor<3> ele(mesh, elm_idx_[i_elm]);
	    arma::mat map_mat = mapping.element_map(ele);
	    for (unsigned int i_point=0; i_point<eval_points->size(); ++i_point) {
	    	global_coords_.get<3, 1>(i_elm*eval_points->size()+i_point)
	            = mapping.project_unit_to_real(RefElement<elemdim>::local_to_bary(eval_points->local_point(i_point)), map_mat);
	    }
    }
}



#endif /* FIELD_VALUE_CACHE_IMPL_HH_ */
