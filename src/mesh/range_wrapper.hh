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
 * @file    range_wrapper.hh
 * @brief   Implementation of range helper class.
 */

#ifndef RANGE_WRAPPER_HH_
#define RANGE_WRAPPER_HH_

#include "mesh/mesh.h"
#include "tools/general_iterator.hh"

/**
 * @brief Range helper class.
 *
 * Allow iterate in bounds given by begin and end range. Class can be used for iterable accessor classes.
 */
template<class Object, class Source>
class Range
{
public:
	Range(const Source * src, unsigned int begin, unsigned int end)
	: source_(src), begin_(begin), end_(end) {
		ASSERT_LE(begin, end).error("Invalid range, begin is greater than end!");
	}

	Iter<Object> begin() {
		return make_iter<Object>( Object(source_, begin_) );
	}

	Iter<Object> end() {
		return make_iter<Object>( Object(source_, end_) );
	}

	inline unsigned int size() const {
		return end_ - begin_;
	}
private:
	const Source * source_;
	unsigned int begin_;
	unsigned int end_;
};

#endif // RANGE_WRAPPER_HH_
