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
 * @brief List of base attributes of Flow123D application.
 *
 * TODO: This structure will be move to flow namespace when it will be available.
 */
struct FlowAttributes {
	/// Actual version of application.
	string version;
	/// Actual revision of application.
	string revision;
	/// Display name of branch.
	string branch;
	/// Display url.
	string url;
};

static FlowAttributes flow_attributes;

namespace Input {
namespace Type {

/**
 * @brief List of common attributes of Input::Type objects.
 *
 * These attributes can be used in any Input::Type object.
 */
struct InputAttributes {
	/**
	 * Indicates that formatter should make the type documentation part
	 * of the documentation of the type that use it. E.g. documentation
	 * of a record key contain documentation of its type.
	 *
	 * Format of value: string
	 */
	string embedded_doc;
	/**
	 * Propose custom target name for hypertext reference.
	 *
	 * Format of value: string
	 */
	string link_name;
	/**
	 * Obsolete type.
	 *
	 * Format of value: bool
	 */
	string obsolete;
	/**
	 * JSON with description of move of the particular type/key (only if
	 * we allow attributes of keys).
	 *
	 * Format of value: string
	 */
	string ist_change;
	/**
	 * Particular for GeoMop, units in machine readable form
	 *
	 * Format of value: string
	 */
	string units;
};

static InputAttributes input_attributes;


} // closing namespace Type
} // closing namespace Input

#endif /* ATTRIBUTE_LIB_HH_ */
