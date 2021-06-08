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

#include "tools/general_iterator.hh"

/**
 * @brief Range helper class.
 *
 * Allow iterate in bounds given by begin and end iterator. Class can be used for iterable accessor classes.
 *
 * Template argument:
 *  - ObjectIn  Type over its instances is iterated,
 *  - ObjectOut Operators '*' and '->' returns objects of this type.
 *
 * Require the template object to implement:
 *  - ObjectIn must be implicitly convertible to ObjectOut type.
 */
template<class ObjectIn, class ObjectOut>
class RangeConvert
{
public:
	/// Constructor.
	RangeConvert(IterConvert<ObjectIn, ObjectOut> begin, IterConvert<ObjectIn, ObjectOut> end)
	: begin_(begin), end_(end) {}

	/// Iterator to begin item of range.
	IterConvert<ObjectIn, ObjectOut> begin() {
		return begin_;
	}

	/// Iterator to end item of range.
	IterConvert<ObjectIn, ObjectOut> end() {
		return end_;
	}

private:
	IterConvert<ObjectIn, ObjectOut> begin_;
	IterConvert<ObjectIn, ObjectOut> end_;
};


/**
 * @brief Range helper class.
 *
 * Same as previous but doesn't provide specialization of operators '*' and '->'.
 */
template<class Object>
class Range
{
public:
	/// Constructor.
	Range(Iter<Object> begin, Iter<Object> end)
	: begin_(begin), end_(end) {}

	/// Iterator to begin item of range.
	Iter<Object> begin() {
		return begin_;
	}

	/// Iterator to end item of range.
	Iter<Object> end() {
		return end_;
	}

private:
	Iter<Object> begin_;
	Iter<Object> end_;
};


#endif // RANGE_WRAPPER_HH_
