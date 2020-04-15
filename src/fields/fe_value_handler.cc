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
template<int rank, int spacedim, class Value>
class FEShapeHandler {
public:

	inline static typename Value::return_type fe_value(FEValues<spacedim> &fe_val, unsigned int i_dof, unsigned int i_qp, unsigned int comp_index)
	{
		ASSERT(false).error("Unsupported format of FieldFE!\n");
		typename Value::return_type ret;
		Value val(ret);
		val.zeros();
		return ret;
	}
};


/// Partial template specialization of FEShapeHandler for scalar fields
template<int spacedim, class Value>
class FEShapeHandler<0, spacedim, Value> {
public:
	inline static typename Value::return_type fe_value(FEValues<3> &fe_val, unsigned int i_dof, unsigned int i_qp, FMT_UNUSED unsigned int comp_index)
	{
		return fe_val.scalar_view(comp_index).value(i_dof, i_qp);
	}
};


/// Partial template specialization of FEShapeHandler for vector fields
template<int spacedim, class Value>
class FEShapeHandler<1, spacedim, Value> {
public:
	inline static typename Value::return_type fe_value(FEValues<3> &fe_val, unsigned int i_dof, unsigned int i_qp, unsigned int comp_index)
	{
		return fe_val.vector_view(comp_index).value(i_dof, i_qp);
	}
};


/// Partial template specialization of FEShapeHandler for tensor fields
template<int spacedim, class Value>
class FEShapeHandler<2, spacedim, Value> {
public:
	inline static typename Value::return_type fe_value(FEValues<3> &fe_val, unsigned int i_dof, unsigned int i_qp, unsigned int comp_index)
	{
		return fe_val.tensor_view(comp_index).value(i_dof, i_qp);
	}
};



template <int elemdim, int spacedim, class Value>
FEValueHandler<elemdim, spacedim, Value>::FEValueHandler()
: value_(r_value_)
{}


template <int elemdim, int spacedim, class Value>
void FEValueHandler<elemdim, spacedim, Value>::initialize(FEValueInitData init_data)
{
	if (dh_)
		WarningOut() << "Multiple initialization of FEValueHandler!";

	dh_ = init_data.dh;
	data_vec_ = init_data.data_vec;
    value_.set_n_comp(init_data.n_comp);
    comp_index_ = init_data.comp_index;
}


template <int elemdim, int spacedim, class Value> inline
typename Value::return_type const &FEValueHandler<elemdim, spacedim, Value>::value(const Point &p, const ElementAccessor<spacedim> &elm)
{
	Armor::array point_list(spacedim, 1, 1);
	point_list.set(0) = Armor::ArmaVec<double,spacedim>( p );
	std::vector<typename Value::return_type> v_list;
	v_list.push_back(r_value_);
	this->value_list(point_list, elm, v_list);
	this->r_value_ = v_list[0];
	return this->r_value_;
}


template <int elemdim, int spacedim, class Value>
void FEValueHandler<elemdim, spacedim, Value>::value_list(const Armor::array  &point_list, const ElementAccessor<spacedim> &elm,
                   std::vector<typename Value::return_type> &value_list)
{
    ASSERT_EQ( point_list.size(), value_list.size() ).error();

	const DHCellAccessor cell = dh_->cell_accessor_from_element( elm.idx() );
	LocDofVec loc_dofs;
	if (boundary_dofs_) loc_dofs = this->get_loc_dof_indices(elm.idx());
	else loc_dofs = cell.get_loc_dof_indices();

	// map points to reference cell, create quadrature and FEValues object
    arma::mat map_mat = MappingP1<elemdim,spacedim>::element_map(elm);
	Quadrature quad(elemdim, point_list.size());
	for (unsigned int k=0; k<point_list.size(); k++)
        quad.set(k) = RefElement<elemdim>::bary_to_local(MappingP1<elemdim,spacedim>::project_real_to_unit(point_list.vec<spacedim>(k), map_mat));
	
	MixedPtr<FiniteElement> fe_mixed_ptr = dh_->ds()->fe(elm);
	std::shared_ptr<FiniteElement<elemdim>> fe_ptr = fe_mixed_ptr.get<elemdim>();
	FEValues<spacedim> fe_values(quad, *fe_ptr, update_values);
    fe_values.reinit( elm );

    for (unsigned int k=0; k<point_list.size(); k++) {
		Value envelope(value_list[k]);
		envelope.zeros();
		for (unsigned int i=0; i<loc_dofs.n_elem; i++) {
			value_list[k] += data_vec_[loc_dofs[i]]
							* FEShapeHandler<Value::rank_, spacedim, Value>::fe_value(fe_values, i, k, comp_index_);
		}
	}
}


template <int elemdim, int spacedim, class Value>
unsigned int FEValueHandler<elemdim, spacedim, Value>::compute_quadrature(std::vector<arma::vec::fixed<3>> & q_points, std::vector<double> & q_weights,
		const ElementAccessor<spacedim> &ele, unsigned int order)
{
	static const double weight_coefs[] = { 1., 1., 2., 6. };

	QGauss qgauss(elemdim, order);
	arma::mat map_mat = MappingP1<elemdim,spacedim>::element_map(ele);

	for(unsigned i=0; i<qgauss.size(); ++i) {
		q_weights[i] = qgauss.weight(i)*weight_coefs[elemdim];
		q_points[i] = MappingP1<elemdim,spacedim>::project_unit_to_real(RefElement<elemdim>::local_to_bary(qgauss.point<elemdim>(i)), map_mat);
	}

	return qgauss.size();
}


template <int spacedim, class Value>
void FEValueHandler<0, spacedim, Value>::initialize(FEValueInitData init_data)
{
	if (dh_)
		WarningOut() << "Multiple initialization of FEValueHandler!";

	dh_ = init_data.dh;
	data_vec_ = init_data.data_vec;
    value_.set_n_comp(init_data.n_comp);
}


template <int spacedim, class Value>
void FEValueHandler<0, spacedim, Value>::value_list(const Armor::array &point_list, const ElementAccessor<spacedim> &elm,
                   std::vector<typename Value::return_type> &value_list)
{
	ASSERT_EQ( point_list.size(), value_list.size() ).error();

	const DHCellAccessor cell = dh_->cell_accessor_from_element( elm.idx() );
	LocDofVec loc_dofs;
	if (boundary_dofs_) loc_dofs = this->get_loc_dof_indices(elm.idx());
	else loc_dofs = cell.get_loc_dof_indices();

	for (unsigned int k=0; k<point_list.size(); k++) {
		Value envelope(value_list[k]);
		envelope.zeros();
		for (unsigned int i=0; i<loc_dofs.n_elem; i++) {
			envelope(i / envelope.n_cols(), i % envelope.n_rows()) += data_vec_[loc_dofs[i]];
		}
	}
}


template <int elemdim, int spacedim, class Value>
FEValueHandler<elemdim, spacedim, Value>::~FEValueHandler()
{
}


// Instantiations of FEValueHandler and FEShapeHandler
#define INSTANCE_VALUE_HANDLER_ALL(dim, spacedim)                                     \
template class FEValueHandler<dim, spacedim, FieldValue<0>::Enum >;                   \
template class FEValueHandler<dim, spacedim, FieldValue<0>::Integer >;                \
template class FEValueHandler<dim, spacedim, FieldValue<0>::Scalar >;                 \
template class FEValueHandler<dim, spacedim, FieldValue<spacedim>::VectorFixed >;     \
template class FEValueHandler<dim, spacedim, FieldValue<spacedim>::TensorFixed >;

#define INSTANCE_VALUE_HANDLER(dim) \
INSTANCE_VALUE_HANDLER_ALL(dim,3)
//INSTANCE_VALUE_HANDLER_ALL(dim,2)

INSTANCE_VALUE_HANDLER(0)
INSTANCE_VALUE_HANDLER(1)
INSTANCE_VALUE_HANDLER(2)
INSTANCE_VALUE_HANDLER(3)

template class FEShapeHandler<0, 3, FieldValue<0>::Enum >;
template class FEShapeHandler<0, 3, FieldValue<0>::Integer >;
template class FEShapeHandler<0, 3, FieldValue<0>::Scalar >;
template class FEShapeHandler<1, 3, FieldValue<3>::VectorFixed >;
template class FEShapeHandler<2, 3, FieldValue<3>::TensorFixed >;
