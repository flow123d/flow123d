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
class Attribute {
public:
    /**
     * This attribute provides a list of names of the free parameters of the subtree of a generic type.
     * Value is a list of strings.
     */
    inline static string generic_parameters()
    { return "_generic_parameters"; }

    /**
     * This attribute is set to value 'true' for the generic types that are the root types of a generic subtree, i.e.
     * The Input::Type::Instance was applied to them.
     * Value is a bool.
     */
    inline static string root_of_generic_subtree()
    { return "_root_of_generic_subtree"; }

    /**
	 * Indicates that formatter should make the type documentation part
	 * of the documentation of the type that use it. E.g. documentation
	 * of a record key contain documentation of its type.
	 *
	 * Value is a bool.
	 */
	inline static string embedded_doc()
	{ return "_embedded_doc"; }

	/**
	 * Propose custom target name for hypertext reference.
	 *
	 * Format of value: string
	 */
	inline static string link_name()
	{ return "_link_name"; }

	/**
	 * Attribute to mark obsolete types. The value of the attribute is the replacement of the
	 * feature, or reason to abandon the type.
	 *
	 * Format of value: bool
	 */
	inline static string obsolete()
	{ return "obsolete"; }

};


} // closing namespace Type
} // closing namespace Input

#endif /* TYPE_ATTRIBUTE_LIB_HH_ */
