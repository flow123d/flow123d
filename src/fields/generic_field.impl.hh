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
	static FE_P<0>::disc fe0(0);
	static FE_P<1>::disc fe1(0);
	static FE_P<2>::disc fe2(0);
	static FE_P<3>::disc fe3(0);
    DOFHandlerMultiDim dh_par(mesh);
    std::shared_ptr<DiscreteSpace> ds = std::make_shared<EqualOrderDiscreteSpace>( &mesh, &fe0, &fe1, &fe2, &fe3);
    dh_par.distribute_dofs(ds);
    std::shared_ptr<DOFHandlerMultiDim> dh = dh_par.sequential();

	auto field_subdomain_data = mesh.get_part()->subdomain_id_field_data();
	std::vector<LongIdx> indices(1);
	VectorMPI *data_vec = new VectorMPI(mesh.n_elements());
	ASSERT_EQ(dh->max_elem_dofs(), 1);
	unsigned int i_ele=0;
	for (auto cell : dh->own_range()) {
		cell.get_loc_dof_indices(indices);
		(*data_vec)[ indices[0] ] = (*field_subdomain_data)[i_ele];
		++i_ele;
	}
    std::shared_ptr< FieldFE<spacedim, DoubleScalar> > field_ptr = std::make_shared< FieldFE<spacedim, DoubleScalar> >();
    field_ptr->set_fe_data(dh);

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
