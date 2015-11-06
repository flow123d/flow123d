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
 * @file    generic_field.hh
 * @brief   
 */

#ifndef GENERIC_FIELDS_HH_
#define GENERIC_FIELDS_HH_

/**
 * @file
 * @brief Fields computed from the mesh data.
 *
 * This file collects fields that are independent of particular equation and
 * depends only on data in mesh.
 */

#include "fields/field.hh"
#include "fields/field_constant.hh"
#include "fields/field_elementwise.hh"

class Mesh;

template <int spacedim>
class GenericField {
public:

	/// Index value type
	typedef typename FieldValue<spacedim>::Integer IntegerScalar;
	/// Index valued field
	typedef Field<spacedim, IntegerScalar> IndexField;

	/**
	 * Returns an instance of a scalar integer field that provides ID's of regions.
	 */
	static auto region_id(Mesh &mesh) -> IndexField;

	/**
	 * Returns an instance of a scalar integer field that provides ID's of subdomains used for
	 * domain decomposition.
	 *
	 * TODO: FieldElementwise just use provided data pointer (unsafe solution), so currently we store the data into a mesh.
	 * Which is safe as long as we have one mesh for whol calculation.
	 * After we have FieldFE that use some sort of Vector class with own memory management, we should use these
	 * to pass the data in safe way.
	 */
	static auto subdomain(Mesh &mesh) -> IndexField;
};


#endif /* GENERAL_FIELDS_HH_ */
