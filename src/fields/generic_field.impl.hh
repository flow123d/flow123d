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
 * @file    generic_field.impl.hh
 * @brief   
 */

#ifndef GENERIC_FIELD_IMPL_HH_
#define GENERIC_FIELD_IMPL_HH_

#include <memory>
#include "mesh/mesh.h"

#include "fields/generic_field.hh"


template <int spacedim>
auto GenericField<spacedim>::region_id(Mesh &mesh) -> IndexField {
	IndexField region_id;
	region_id.name("region_id");
	region_id.units( UnitSI::dimensionless() );
	region_id.set_mesh(mesh);

	RegionSet all_regions=mesh.region_db().get_region_set("ALL");
	for(Region reg : all_regions) {
		auto field_algo=std::make_shared<FieldConstant<spacedim, IntegerScalar>>();
		field_algo->set_value(reg.id());
		region_id.set_field(
				{reg} ,
				field_algo);
	}

	//region_id.set_time(0.0);

	return region_id;
}

template <int spacedim>
auto GenericField<spacedim>::subdomain(Mesh &mesh) -> IndexField {
	auto field_subdomain_data= mesh.get_part()->subdomain_id_field_data();

	IndexField subdomain;
	subdomain.name("subdomain");
	subdomain.units( UnitSI::dimensionless() );
	subdomain.set_mesh(mesh);

	subdomain.set_field(
		mesh.region_db().get_region_set("ALL"),
		make_shared< FieldElementwise<spacedim, FieldValue<3>::Integer> >(field_subdomain_data, 1),
		0.0); // time=0.0


	//subdomain.set_time(0.0);
	return subdomain;
}



#endif /* GENERAL_FIELD_IMPL_HH_ */
