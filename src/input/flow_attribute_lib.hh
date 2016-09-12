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
class FlowAttribute {
public:

	/**
	 * Specify a unit of a field.
	 *
	 * Format of value: record with base units as keys and exponents as their values, e.g.
	 * { "m" : 0, "md" : 0, "kg" : 0, "s" : 0, "A" : 0, "K" : 0, "mol" : 0, "cd" : 0 }
	 */
	inline static string field_unit()
	{ return "key_field_unit"; }

    /**
     * Specify a value of the field.
     *
     * Format of value: record e.g.
     * { "subfields": true, "shape": [ 1, 1 ], "type": "Float" }
     */
    inline static string field_value_shape()
    { return "field_value_shape"; }

    /**
     * Specify path where multifields of the equation get their names.
     *
     * Format of value: string, pattern of address, e.g.
     * "/problem/solute_equation/substances/ * /name"
     */
    inline static string subfields_address()
    { return "subfields_address"; }

};

#endif /* ATTRIBUTE_LIB_HH_ */
