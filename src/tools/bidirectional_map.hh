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
 * @file    bidirectional_map.hh
 * @brief   Implementation of bidirectional map.
 */

#ifndef BIDIRECTIONAL_MAP_HH_
#define BIDIRECTIONAL_MAP_HH_

/** @brief Bidirectional map templated by <T, unsigned int>.
 * Provides bidirectional access between pair of values.
 */
template<typename T>
class BidirectionalMap
{
public:
	BidirectionalMap(unsigned int init_size = 0);

	unsigned int size();

	void set_item(T val, unsigned int idx);

	unsigned int add_item(T val);

	int get_index(T val);

	T operator[](unsigned int idx);

private:
    std::vector<T> vals_vec_;             ///< Space to save values.
    std::map<T, unsigned int> vals_map_;  ///< Maps values to indexes into vals_vec_.
};

// --------------------------------------------------- BidirectionalMap INLINE implementation -----------
template<typename T>
inline BidirectionalMap<T>::BidirectionalMap(unsigned int init_size)
: vals_vec_(init_size, -1)
{}

template<typename T>
inline unsigned int BidirectionalMap<T>::size() {
	return vals_map_.size();
}

template<typename T>
inline void BidirectionalMap<T>::set_item(T val, unsigned int idx) {
	ASSERT_LT( idx, vals_vec_.size() )(idx)(vals_vec_.size()).error("Value id is out of vector size.");
	ASSERT( vals_vec_[idx] != -1 )(idx).error("Repeated setting of item.");
	vals_map_[val] = idx;
	vals_vec_[idx] = val;
}

template<typename T>
inline unsigned int BidirectionalMap<T>::add_item(T val) {
	ASSERT( vals_map_.find(val) == vals_map_.end() )(val).error("Can not add item since it already exists.");
	vals_map_[val] = vals_vec_.size();
	vals_vec_.push_back(val);
	return vals_map_[val];
}

template<typename T>
inline int BidirectionalMap<T>::get_index(T val) {
	typename std::map<T, unsigned int>::iterator iter = vals_map_.find(val);
	if (iter == vals_map_.end()) return -1;
	else return iter.second;
}

template<typename T>
inline T BidirectionalMap<T>::operator[](unsigned int idx) {
	ASSERT( idx < vals_vec_.size() );
	return vals_vec_[idx];
}


#endif // BIDIRECTIONAL_MAP_HH_
