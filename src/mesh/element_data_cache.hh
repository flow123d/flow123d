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
#include "system/tokenizer.hh"


class ElementDataCacheBase {
public:
	/// Constructor.
	ElementDataCacheBase()
	: time_(-std::numeric_limits<double>::infinity()),
	  quantity_name_("") {}

	/// Destructor
	virtual ~ElementDataCacheBase() {}

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

	/**
	 * Read ascii data of given \p i_row from tokenizer
	 */
	virtual void read_ascii_data(Tokenizer &tok, unsigned int n_components, unsigned int i_row)=0;

	/**
	 * Read binary data of given \p i_row from data stream
	 */
	virtual void read_binary_data(std::istream &data_stream, unsigned int n_components, unsigned int i_row)=0;

protected:
	/// time step stored in cache
	double time_;
	/// name of quantity stored in cache
	std::string quantity_name_;
};



struct MeshDataHeader;


template <typename T>
class ElementDataCache : public ElementDataCacheBase {
public:
	typedef std::shared_ptr< std::vector<T> > ComponentDataPtr;
	typedef std::vector< ComponentDataPtr > CacheData;

	/// Default constructor
	ElementDataCache();

	/**
	 * Constructor.
	 *
	 * @param data_header   Set data members time_ and quantity_name_
	 * @param size_of_cache Count of columns of data cache
	 * @param row_vec_size  Count of rows of data cache
	 */
	ElementDataCache(MeshDataHeader data_header, unsigned int size_of_cache, unsigned int row_vec_size);

	/// Return vector of element data for get component.
	ComponentDataPtr get_component_data(unsigned int component_idx);

	/**
	 * Create data cache with given count of columns (\p size_of_cache) and rows (\p row_vec_size).
	 */
	static CacheData create_data_cache(unsigned int size_of_cache, unsigned int row_vec_size);

	/// Implements @p ElementDataCacheBase::read_ascii_data.
	void read_ascii_data(Tokenizer &tok, unsigned int n_components, unsigned int i_row) override;

	/// Implements @p ElementDataCacheBase::read_binary_data.
	void read_binary_data(std::istream &data_stream, unsigned int n_components, unsigned int i_row) override;

protected:
	/**
	 * Table of element data.
	 *
	 * For every components contains vector of element data.
	 */
	CacheData data_;
    
};

#endif /* ELEMENT_DATA_CACHE_HH_ */
