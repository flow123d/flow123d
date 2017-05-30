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
 * @file    element_data_cache.cc
 * @brief
 */


#include "mesh/element_data_cache.hh"
#include "mesh/msh_basereader.hh"
#include "boost/lexical_cast.hpp"



template <typename T>
ElementDataCache<T>::ElementDataCache()
: ElementDataCacheBase() {}


template <typename T>
ElementDataCache<T>::ElementDataCache(MeshDataHeader data_header, unsigned int size_of_cache, unsigned int row_vec_size) {
	this->time_ = data_header.time;
	this->quantity_name_ = data_header.field_name;
	this->data_ = create_data_cache(size_of_cache, row_vec_size);
}


template <typename T>
typename ElementDataCache<T>::ComponentDataPtr ElementDataCache<T>::get_component_data(unsigned int component_idx) {
	ASSERT_LT(component_idx, data_.size()).error("Index of component is out of range.\n");
	return data_[component_idx];
}


template <typename T>
typename ElementDataCache<T>::CacheData ElementDataCache<T>::create_data_cache(unsigned int size_of_cache, unsigned int row_vec_size) {
    typename ElementDataCache<T>::CacheData data_cache(size_of_cache);
    for (unsigned int i=0; i<size_of_cache; ++i) {
		typename ElementDataCache<T>::ComponentDataPtr row_vec = std::make_shared<std::vector<T>>();
		row_vec->resize(row_vec_size);
		data_cache[i] = row_vec;
    }

    return data_cache;
}


template <typename T>
void ElementDataCache<T>::read_ascii_data(Tokenizer &tok, unsigned int n_components, unsigned int i_row) {
	unsigned int idx;
	for (unsigned int i_vec=0; i_vec<data_.size(); ++i_vec) {
		idx = i_row * n_components;
		std::vector<T> &vec = *( data_[i_vec].get() );
		for (unsigned int i_col=0; i_col < n_components; ++i_col, ++idx) {
			vec[idx] = boost::lexical_cast<T>(*tok);
			++tok;
		}
	}
}


template <typename T>
void ElementDataCache<T>::read_binary_data(std::istream &data_stream, unsigned int n_components, unsigned int i_row) {
	unsigned int idx;
	for (unsigned int i_vec=0; i_vec<data_.size(); ++i_vec) {
		idx = i_row * n_components;
		std::vector<T> &vec = *( data_[i_vec].get() );
		for (unsigned int i_col=0; i_col < n_components; ++i_col, ++idx) {
			data_stream.read(reinterpret_cast<char *>(&vec[idx]), sizeof(T));
		}
	}
}


// explicit instantiation of template class
template class ElementDataCache<unsigned int>;
template class ElementDataCache<int>;
template class ElementDataCache<double>;
