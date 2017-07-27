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
#include "mesh/bounding_box.hh"


/**
 * Helper class, allow to simplify computing value of FieldFE.
 *
 * Use correct method FEValues<...>::shape_xxx given with Value::rank_.
 * Practical use have only instances with rank template parameters 0 and 1 (Scalar and Vector Fields, see below).
 */
template<int rank, int elemdim, int spacedim, class Value>
class FEShapeHandler {
public:
	inline static typename Value::return_type fe_value(FEValues<elemdim,3> &fe_val, unsigned int i_dof, unsigned int i_qp)
	{
		ASSERT(false).error("Unsupported format of FieldFE!\n");
		typename Value::return_type ret;
		Value val(ret);
		val.zeros();
		return ret;
	}
};


/// Partial template specialization of FEShapeHandler for scalar fields
template<int elemdim, int spacedim, class Value>
class FEShapeHandler<0, elemdim, spacedim, Value> {
public:
	inline static typename Value::return_type fe_value(FEValues<elemdim,3> &fe_val, unsigned int i_dof, unsigned int i_qp)
	{
		return fe_val.shape_value(i_dof, i_qp);
	}
};


/// Partial template specialization of FEShapeHandler for vector fields
template<int elemdim, int spacedim, class Value>
class FEShapeHandler<1, elemdim, spacedim, Value> {
public:
	inline static typename Value::return_type fe_value(FEValues<elemdim,3> &fe_val, unsigned int i_dof, unsigned int i_qp)
	{
		return fe_val.shape_vector(i_dof, i_qp);
	}
};



template <int elemdim, int spacedim, class Value>
FEValueHandler<elemdim, spacedim, Value>::FEValueHandler()
: value_(r_value_),
  map_(nullptr)
{}


template <int elemdim, int spacedim, class Value>
void FEValueHandler<elemdim, spacedim, Value>::initialize(FEValueInitData init_data, MappingP1<elemdim,3> *map)
{
	ASSERT_EQ(dof_indices.size(), 0).error("Multiple initialization.");

	dh_ = init_data.dh;
	data_vec_ = init_data.data_vec;
    dof_indices.resize(init_data.ndofs);
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

    arma::mat map_mat = map_->element_map(*elm.element());
	for (unsigned int k=0; k<point_list.size(); k++) {
		Quadrature<elemdim> quad(1);
        quad.set_point(0, RefElement<elemdim>::bary_to_local(map_->project_real_to_unit(point_list[k], map_mat)));

        FEValues<elemdim,3> fe_values(*map_, quad, *dh_->fe<elemdim>(), update_values);
		fe_values.reinit(cell);

		Value envelope(value_list[k]);
		envelope.zeros();
		for (unsigned int i=0; i<dh_->fe<elemdim>()->n_dofs(); i++)
			value_list[k] += (*data_vec_)[dof_indices[i]]
										  * FEShapeHandler<Value::rank_, elemdim, spacedim, Value>::fe_value(fe_values, i, 0);
	}
}


template <int elemdim, int spacedim, class Value>
bool FEValueHandler<elemdim, spacedim, Value>::contains_point(arma::vec point, Element &elm)
{
	ASSERT_PTR(map_).error();

	arma::vec projection = map_->project_real_to_unit(point, map_->element_map(elm));
	return (projection.min() >= -BoundingBox::epsilon);
}


template <int elemdim, int spacedim, class Value>
FEValueHandler<elemdim, spacedim, Value>::~FEValueHandler()
{}


// Instantiations of FEValueHandler and FEShapeHandler
#define INSTANCE_VALUE_HANDLER_ALL(dim, spacedim)                                     \
template class FEValueHandler<dim, spacedim, FieldValue<0>::Enum >;                   \
template class FEValueHandler<dim, spacedim, FieldValue<0>::Integer >;                \
template class FEValueHandler<dim, spacedim, FieldValue<0>::Scalar >;                 \
template class FEValueHandler<dim, spacedim, FieldValue<spacedim>::VectorFixed >;     \
template class FEValueHandler<dim, spacedim, FieldValue<spacedim>::TensorFixed >;     \
template class FEShapeHandler<0, dim, spacedim, FieldValue<0>::Enum >;                \
template class FEShapeHandler<0, dim, spacedim, FieldValue<0>::Integer >;             \
template class FEShapeHandler<0, dim, spacedim, FieldValue<0>::Scalar >;              \
template class FEShapeHandler<1, dim, spacedim, FieldValue<spacedim>::VectorFixed >;  \
template class FEShapeHandler<2, dim, spacedim, FieldValue<spacedim>::TensorFixed >;

#define INSTANCE_VALUE_HANDLER(dim) \
INSTANCE_VALUE_HANDLER_ALL(dim,2)   \
INSTANCE_VALUE_HANDLER_ALL(dim,3)

INSTANCE_VALUE_HANDLER(1);
INSTANCE_VALUE_HANDLER(2);
INSTANCE_VALUE_HANDLER(3);
