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
#include "mesh/partitioning.hh"
#include "mesh/accessors.hh"

#include "fields/generic_field.hh"
#include "fields/field_fe.hh"
#include "la/vector_mpi.hh"
#include "fields/fe_value_handler.hh"

#include "fem/mapping_p1.hh"
#include "fem/fe_p.hh"
#include "fem/dofhandler.hh"
#include "fem/discrete_space.hh"


template <int spacedim>
auto GenericField<spacedim>::region_id(Mesh &mesh) -> IndexField {
	IndexField region_id;
	region_id.name("region_id");
	region_id.units( UnitSI::dimensionless() );
	region_id.set_mesh(mesh);

	RegionSet all_regions=mesh.region_db().get_region_set("ALL");
	for(Region reg : all_regions) {
		auto field_algo=std::make_shared<FieldConstant<spacedim, DoubleScalar>>();
		field_algo->set_value(reg.id());
		region_id.set_field(
				{reg} ,
				field_algo,
				0.0); // time=0.0
	}
	return region_id;
}

template <int spacedim>
auto GenericField<spacedim>::subdomain(Mesh &mesh) -> IndexField {
    MixedPtr<FE_P_disc> fe(0);
    std::shared_ptr<DOFHandlerMultiDim> dh = std::make_shared<DOFHandlerMultiDim>(mesh);
    std::shared_ptr<DiscreteSpace> ds = std::make_shared<EqualOrderDiscreteSpace>( &mesh, fe);
    dh->distribute_dofs(ds);

	auto field_subdomain_data = mesh.get_part()->subdomain_id_field_data();
	unsigned int data_size = field_subdomain_data->size();
	VectorMPI data_vec(data_size);
	ASSERT_EQ(dh->max_elem_dofs(), 1);
	unsigned int i_ele=0;
	for (auto cell : dh->own_range()) {
		data_vec[ cell.get_loc_dof_indices()(0) ] = (*field_subdomain_data)[i_ele];
		++i_ele;
	}
    std::shared_ptr< FieldFE<spacedim, DoubleScalar> > field_ptr = std::make_shared< FieldFE<spacedim, DoubleScalar> >();
    field_ptr->set_fe_data(dh, 0, data_vec);

	IndexField subdomain;
	subdomain.name("subdomain");
	subdomain.units( UnitSI::dimensionless() );
	subdomain.set_mesh(mesh);

    subdomain.set_field(
		mesh.region_db().get_region_set("ALL"),
		field_ptr,
		0.0); // time=0.0

	return subdomain;
}



#endif /* GENERAL_FIELD_IMPL_HH_ */
