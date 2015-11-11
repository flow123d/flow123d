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
		ASSERT(component_idx < data_.size(), "Index of component is out of range.\n");
		return data_[component_idx];
	}

protected:
	/**
	 * Table of element data.
	 *
	 * For every components contains vector of element data.
	 */
	CacheData data_;
    
};

#endif /* ELEMENT_DATA_CACHE_HH_ */
