/*
 * generic_field.impl.hh
 *
 *  Created on: Dec 9, 2014
 *      Author: jb
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
