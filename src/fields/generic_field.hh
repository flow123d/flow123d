/*
 * general_fields.hh
 *
 *  Created on: Dec 9, 2014
 *      Author: jb
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
