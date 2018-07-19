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
 * @file    generic_interpolator.cc
 * @brief
 */

#include "fields/generic_interpolator.hh"
#include "fields/field_instances.hh"	// for instantiation macros
#include "fields/field_fe.hh"
#include "io/reader_cache.hh"


template <int spacedim, class Value>
GenericInterpolator<spacedim, Value>::GenericInterpolator()
{}

template <int spacedim, class Value>
void GenericInterpolator<spacedim, Value>::interpolate(FieldFE<spacedim, Value> & field_out, FieldFE<spacedim, Value> & field_in)
{
	std::shared_ptr<Mesh> source_mesh = ReaderCache::get_mesh(field_in.reader_file_);
	std::vector<double> sum_val(4);
	std::vector<unsigned int> elem_count(4);
	std::vector<unsigned int> searched_elements; // stored suspect elements in calculating the intersection

	for (auto ele : field_out.dh_->mesh()->elements_range()) {
		searched_elements.clear();
		source_mesh->get_bih_tree().find_point(ele.centre(), searched_elements);
		std::fill(sum_val.begin(), sum_val.end(), 0.0);
		std::fill(elem_count.begin(), elem_count.end(), 0);
		for (std::vector<unsigned int>::iterator it = searched_elements.begin(); it!=searched_elements.end(); it++) {
			ElementAccessor<3> elm = source_mesh->element_accessor(*it);
			bool contains=false;
			switch (elm->dim()) {
			case 1:
				contains = field_in.value_handler1_.get_mapping()->contains_point(ele.centre(), elm);
				break;
			case 2:
				contains = field_in.value_handler2_.get_mapping()->contains_point(ele.centre(), elm);
				break;
			case 3:
				contains = field_in.value_handler3_.get_mapping()->contains_point(ele.centre(), elm);
				break;
			default:
				ASSERT(false).error("Invalid element dimension!");
			}
			if (contains) {
				// projection point in element
				sum_val[elm->dim()] += (*field_in.data_vec_)[*it];
				++elem_count[elm->dim()];
			}
		}
		unsigned int dim = ele->dim(); // dim+1 for boundary
		double elem_value = 0.0;
		do {
			if (elem_count[dim] > 0) {
				elem_value = sum_val[dim] / elem_count[dim];
				break;
			}
			++dim;
		} while (dim<4);

		std::cout << " - " << field_out.dh_->get_dof_indices( ele, field_out.dof_indices_);
		ASSERT_LT_DBG( field_out.dof_indices_[0], (int)field_out.data_vec_->size());
		(*field_out.data_vec_)[field_out.dof_indices_[0]] = elem_value * field_in.unit_conversion_coefficient_;
	}
	std::cout << std::endl;
}

// Instantiations of GenericInterpolator
INSTANCE_ALL(GenericInterpolator)
