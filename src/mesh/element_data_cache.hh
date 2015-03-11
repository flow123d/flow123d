/*
 * element_data_cache.hh
 *
 *  Created on: Jan 28, 2013
 *      Author: jb
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
