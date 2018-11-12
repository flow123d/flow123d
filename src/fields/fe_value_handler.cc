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
#include "la/vector_mpi.hh"
#include "fem/mapping_p1.hh"
#include "fem/fe_values.hh"
#include "quadrature/quadrature.hh"
#include "quadrature/quadrature_lib.hh"
#include "mesh/bounding_box.hh"
#include "mesh/accessors.hh"
#include "fem/fe_values_views.hh"
#include "fem/dh_cell_accessor.hh"


/**
 * Helper class, allow to simplify computing value of FieldFE.
 *
 * Use correct method FEValues<...>::shape_xxx given with Value::rank_.
 * Is done by class partial specialization as, we were not able to do this using function overloading (since
 * they differ only by return value) and partial specialization of the function templates is not supported  in C++.
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
		return fe_val.vector_view(0).value(i_dof, i_qp);
	}
};


/// Partial template specialization of FEShapeHandler for tensor fields
template<int elemdim, int spacedim, class Value>
class FEShapeHandler<2, elemdim, spacedim, Value> {
public:
	inline static typename Value::return_type fe_value(FEValues<elemdim,3> &fe_val, unsigned int i_dof, unsigned int i_qp)
	{
		return fe_val.tensor_view(0).value(i_dof, i_qp);
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
	if (dof_indices.size() > 0)
		WarningOut() << "Multiple initialization of FEValueHandler!";

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

    ElementAccessor<3> cell = dh_->mesh()->element_accessor( elm.mesh_idx() ); // non-const ElementAccessor
    if (boundary_dofs_) this->get_dof_indices( elm, dof_indices);
    else dh_->get_dof_indices( cell, dof_indices );

    arma::mat map_mat = map_->element_map(elm);
    for (unsigned int k=0; k<point_list.size(); k++) {
		Quadrature<elemdim> quad(1);
        quad.set_point(0, RefElement<elemdim>::bary_to_local(map_->project_real_to_unit(point_list[k], map_mat)));

		FEValues<elemdim,3> fe_values(*this->get_mapping(), quad, *dh_->fe<elemdim>(cell), update_values);
		fe_values.reinit( const_cast<ElementAccessor<spacedim> &>(elm) );

		Value envelope(value_list[k]);
		envelope.zeros();
		for (unsigned int i=0; i<dh_->fe<elemdim>(cell)->n_dofs(); i++) {
			value_list[k] += (*data_vec_)[dof_indices[i]]
										  * FEShapeHandler<Value::rank_, elemdim, spacedim, Value>::fe_value(fe_values, i, 0);
		}
	}
}


template <int elemdim, int spacedim, class Value>
unsigned int FEValueHandler<elemdim, spacedim, Value>::compute_quadrature(std::vector<arma::vec::fixed<3>> & q_points, std::vector<double> & q_weights,
		const ElementAccessor<spacedim> &ele, unsigned int order)
{
	static const double weight_coefs[] = { 1., 1., 2., 6. };

	QGauss<elemdim> qgauss(order);
	arma::mat map_mat = map_->element_map(ele);

	for(unsigned i=0; i<qgauss.size(); ++i) {
		q_weights[i] = qgauss.weight(i)*weight_coefs[elemdim];
		q_points[i] = map_->project_unit_to_real(RefElement<elemdim>::local_to_bary(qgauss.point(i)), map_mat);
	}

	return qgauss.size();
}


template <int elemdim, int spacedim, class Value>
unsigned int FEValueHandler<elemdim, spacedim, Value>::get_dof_indices(const ElementAccessor<3> &cell, std::vector<LongIdx> &indices) const
{
    unsigned int ndofs = this->value_.n_rows() * this->value_.n_cols();
    for (unsigned int k=0; k<ndofs; k++) {
        indices[k] = (*boundary_dofs_)[ndofs*cell.idx()+k];
    }
    return ndofs;
}


template <int spacedim, class Value>
void FEValueHandler<0, spacedim, Value>::initialize(FEValueInitData init_data)
{
	if (dof_indices.size() > 0)
		WarningOut() << "Multiple initialization of FEValueHandler!";

	dh_ = init_data.dh;
	data_vec_ = init_data.data_vec;
    dof_indices.resize(init_data.ndofs);
    value_.set_n_comp(init_data.n_comp);
}


template <int spacedim, class Value>
void FEValueHandler<0, spacedim, Value>::value_list(const std::vector< Point >  &point_list, const ElementAccessor<spacedim> &elm,
                   std::vector<typename Value::return_type> &value_list)
{
	ASSERT_EQ( point_list.size(), value_list.size() ).error();

	ElementAccessor<3> cell = dh_->mesh()->element_accessor( elm.mesh_idx() ); // non-const ElementAccessor
	if (boundary_dofs_) this->get_dof_indices( elm, dof_indices);
	else dh_->get_dof_indices( cell, dof_indices );

	for (unsigned int k=0; k<point_list.size(); k++) {
		Value envelope(value_list[k]);
		envelope.zeros();
		for (unsigned int i=0; i<dh_->fe<0>(cell)->n_dofs(); i++) {
			envelope(i / envelope.n_cols(), i % envelope.n_rows()) += (*data_vec_)[dof_indices[i]];
		}
	}
}


template <int spacedim, class Value>
unsigned int FEValueHandler<0, spacedim, Value>::get_dof_indices(const ElementAccessor<3> &cell, std::vector<LongIdx> &indices) const
{
    unsigned int ndofs = this->value_.n_rows() * this->value_.n_cols();
    for (unsigned int k=0; k<ndofs; k++) {
        indices[k] = (*boundary_dofs_)[ndofs*cell.idx()+k];
    }
    return ndofs;
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
INSTANCE_VALUE_HANDLER_ALL(dim,3)
//INSTANCE_VALUE_HANDLER_ALL(dim,2)   \

INSTANCE_VALUE_HANDLER(0);
INSTANCE_VALUE_HANDLER(1);
INSTANCE_VALUE_HANDLER(2);
INSTANCE_VALUE_HANDLER(3);
