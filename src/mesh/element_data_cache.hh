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
 * @file    element_data_cache.hh
 * @brief   
 */

#ifndef ELEMENT_DATA_CACHE_HH_
#define ELEMENT_DATA_CACHE_HH_

#include <vector>
#include <string>
#include <memory>
#include "system/system.hh"


class ElementDataCacheBase {
public:
	/// Constructor.
	ElementDataCacheBase()
	: time_(-std::numeric_limits<double>::infinity()),
	  quantity_name_("") {}

	/// Getter for time of cache
	double get_time()
	{ return time_; }

	/// Getter for quantity name of cache
	std::string get_quantity_name()
	{ return quantity_name_; }

	/// Check if cache stored actual data
	bool is_actual(double time, std::string quantity_name) {
		return (time_ == time) && (quantity_name_ == quantity_name);
	}

protected:
	/// time step stored in cache
	double time_;
	/// name of quantity stored in cache
	std::string quantity_name_;
};


template <typename T>
class ElementDataCache : public ElementDataCacheBase {
public:
	typedef std::shared_ptr< std::vector<T> > ComponentDataPtr;
	typedef std::vector< ComponentDataPtr > CacheData;

	/// Constructor.
	ElementDataCache(double time, std::string quantity_name, CacheData data)
	: data_(data) {
		this->time_ = time;
		this->quantity_name_ = quantity_name;
	}

	/// Return vector of element data for get component.
	ComponentDataPtr get_component_data(unsigned int component_idx) {
		ASSERT_LT(component_idx, data_.size()).error("Index of component is out of range.\n");
		return data_[component_idx];
	}

	static CacheData create_data_cache(unsigned int size_of_cache, unsigned int row_vec_size) {
	    typename ElementDataCache<T>::CacheData data_cache(size_of_cache);
	    for (unsigned int i=0; i<size_of_cache; ++i) {
			typename ElementDataCache<T>::ComponentDataPtr row_vec = std::make_shared<std::vector<T>>();
			row_vec->resize(row_vec_size);
			data_cache[i] = row_vec;
	    }

	    return data_cache;
	}

protected:
	/// Empty constructor accessible only for descendants.
	ElementDataCache()
	{}
	/**
	 * Table of element data.
	 *
	 * For every components contains vector of element data.
	 */
	CacheData data_;
    
};

#endif /* ELEMENT_DATA_CACHE_HH_ */
