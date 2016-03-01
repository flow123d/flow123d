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
 * @file    attribute_lib.hh
 * @brief
 */

#ifndef ATTRIBUTE_LIB_HH_
#define ATTRIBUTE_LIB_HH_

using namespace std;

/**
 * @brief Class with static methods provided special attributes of Flow123D application.
 *
 * This class contains only attributes typical for Flow123D application. Base attributes
 * are stored in @p Input::Type::Attributes.
 */
class FlowAttributes {
	/**
	 * Reference to generic type from which Input::Type object is derived.
	 *
	 * Format of value: hash of generic type
	 */
	inline static string generic_type()
	{ return "generic_type"; }
	/**
	 * List of parameters used in generic types or their instances
	 *
	 * Format of value: list of names or list of pairs (name : value)
	 */
	inline static string parameters()
	{ return "parameters"; }
	/**
	 * Particular for GeoMop, units in machine readable form
	 *
	 * Format of value: string
	 */
	inline static string units()
	{ return "units"; }
};

#endif /* ATTRIBUTE_LIB_HH_ */
