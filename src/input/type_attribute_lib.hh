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
 * @file    type_attribute_lib.hh
 * @brief
 */

#ifndef TYPE_ATTRIBUTE_LIB_HH_
#define TYPE_ATTRIBUTE_LIB_HH_

using namespace std;

namespace Input {
namespace Type {

/**
 * @brief Class with static methods provided common attributes of Input::Type objects.
 *
 * These attributes can be used in any Input::Type object.
 */
class Attributes {
public:
	/**
	 * Indicates that formatter should make the type documentation part
	 * of the documentation of the type that use it. E.g. documentation
	 * of a record key contain documentation of its type.
	 *
	 * Format of value: string
	 */
	inline static string embedded_doc()
	{ return "embedded_doc"; }
	/**
	 * Propose custom target name for hypertext reference.
	 *
	 * Format of value: string
	 */
	inline static string link_name()
	{ return "link_name"; }
	/**
	 * Obsolete type.
	 *
	 * Format of value: bool
	 */
	inline static string obsolete()
	{ return "obsolete"; }
	/**
	 * JSON with description of move of the particular type/key (only if
	 * we allow attributes of keys).
	 *
	 * Format of value: string
	 */
	inline static string ist_change()
	{ return "ist_change"; }
};


} // closing namespace Type
} // closing namespace Input

#endif /* TYPE_ATTRIBUTE_LIB_HH_ */
