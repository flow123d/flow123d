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
 * @file    fe_value_handler.cc
 * @brief
 */

#include "fields/fe_value_handler.hh"
#include "fem/mapping_p1.hh"
#include "fem/fe_values.hh"
#include "quadrature/quadrature.hh"


template <int elemdim, int spacedim, class Value>
FEValueHandler<elemdim, spacedim, Value>::FEValueHandler()
: dof_indices(nullptr),
  value_(r_value_),
  map_(nullptr)
{}


template <int elemdim, int spacedim, class Value>
void FEValueHandler<elemdim, spacedim, Value>::initialize(FEValueInitData init_data, Mapping<elemdim,3> *map)
{
	ASSERT(dof_indices == nullptr).error("Multiple initialization.");

	dh_ = init_data.dh;
	data_vec_ = init_data.data_vec;
    dof_indices = new unsigned int[init_data.ndofs];
    value_.set_n_comp(init_data.n_comp);

    if (map == nullptr) {
		// temporary solution - these objects will be set through FieldCommon
		map_ = new MappingP1<elemdim,3>();
    } else {
    	map_ = map;
    }
}


template <int elemdim, int spacedim, class Value> inline
typename Value::return_type const &FEValueHandler<elemdim, spacedim, Value>::value(const Point &p, const ElementAccessor<spacedim> &elm)
{
	std::vector<Point> point_list;
	point_list.push_back(p);
	std::vector<typename Value::return_type> v_list;
	v_list.push_back(r_value_);
	this->value_list(point_list, elm, v_list);
	this->r_value_ = v_list[0];
	return this->r_value_;
}


template <int elemdim, int spacedim, class Value>
void FEValueHandler<elemdim, spacedim, Value>::value_list(const std::vector< Point >  &point_list, const ElementAccessor<spacedim> &elm,
                   std::vector<typename Value::return_type> &value_list)
{
	ASSERT_PTR(map_).error();
	ASSERT_EQ( point_list.size(), value_list.size() ).error();

    DOFHandlerBase::CellIterator cell = dh_->mesh()->element( elm.idx() );
	dh_->get_loc_dof_indices(cell, dof_indices);

	arma::mat::fixed<3,elemdim> m;
	for (unsigned i=0; i<elemdim; ++i) {
		m.col(i) = elm.element()->node[i+1]->point() - elm.element()->node[0]->point();
	}
	arma::mat::fixed<elemdim,3> im = pinv(m);

	for (unsigned int k=0; k<point_list.size(); k++) {
		Point p_rel = point_list[k] - elm.element()->node[0]->point();
		Quadrature<elemdim> quad(1);
		quad.set_point(0, im*p_rel);

		FEValues<elemdim,3> fe_values(*this->get_mapping(), quad, *dh_->fe<elemdim>(), update_values);
		fe_values.reinit(cell);

		Value envelope(value_list[k]);
		envelope.zeros();

		if (dh_->fe<elemdim>()->is_scalar()) {
			for (unsigned int i=0; i<dh_->fe<elemdim>()->n_dofs(); i++)
				value_list[k] += (*data_vec_)[dof_indices[i]]*fe_values.shape_value(i, 0);
		}
		else {
			arma::vec3 value;
			value.zeros();
			for (unsigned int i=0; i<dh_->fe<elemdim>()->n_dofs(); i++)
				value += (*data_vec_)[dof_indices[i]]*fe_values.shape_vector(i, 0);
			for (int i=0; i<3; i++)
				envelope(i,0) = value(i);
		}
	}
}


template <int elemdim, int spacedim, class Value>
FEValueHandler<elemdim, spacedim, Value>::~FEValueHandler()
{
	if (dof_indices != nullptr) delete[] dof_indices;
}


// Instantiations of FEValueHandler
#define INSTANCE_VALUE_HANDLER_ALL(dim, spacedim)                            \
template class FEValueHandler<dim, spacedim, FieldValue<0>::Enum >;          \
template class FEValueHandler<dim, spacedim, FieldValue<0>::Integer >;       \
template class FEValueHandler<dim, spacedim, FieldValue<0>::Scalar >;        \
template class FEValueHandler<dim, spacedim, FieldValue<2>::VectorFixed >;   \
template class FEValueHandler<dim, spacedim, FieldValue<2>::TensorFixed >;   \
template class FEValueHandler<dim, spacedim, FieldValue<3>::VectorFixed >;   \
template class FEValueHandler<dim, spacedim, FieldValue<3>::TensorFixed >;

#define INSTANCE_VALUE_HANDLER(dim) \
INSTANCE_VALUE_HANDLER_ALL(dim,2)   \
INSTANCE_VALUE_HANDLER_ALL(dim,3)

INSTANCE_VALUE_HANDLER(1);
INSTANCE_VALUE_HANDLER(2);
INSTANCE_VALUE_HANDLER(3);
