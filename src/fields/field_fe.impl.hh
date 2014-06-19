/*
 * field_constant.impl.hh
 *
 *  Created on: Dec 15, 2012
 *      Author: jb
 */

#ifndef FIELD_FE_IMPL_HH_
#define FIELD_FE_IMPL_HH_

#include "fields/field_fe.hh"
#include "input/input_type.hh"
#include "quadrature/quadrature.hh"
#include "fem/fe_values.hh"
#include "fem/finite_element.hh"




/// Implementation.

namespace it = Input::Type;




template <int spacedim, class Value>
const int FieldFE<spacedim, Value>::registrar =
		Input::Factory<FactoryBaseType, unsigned int>::template register_class< FieldFE<spacedim, Value> >("FieldFE");



template <int spacedim, class Value>
FieldFE<spacedim, Value>::FieldFE( unsigned int n_comp)
: FieldAlgorithmBase<spacedim, Value>(n_comp),
  dh_(nullptr),
  data_(nullptr),
  data_vec_(nullptr),
  dof_indices(nullptr),
  map1_(nullptr),
  map2_(nullptr),
  map3_(nullptr)
{}


template <int spacedim, class Value>
void FieldFE<spacedim, Value>::set_fe_data(const DOFHandlerMultiDim *dh,
		Mapping<1,3> *map1,
		Mapping<2,3> *map2,
		Mapping<3,3> *map3,
		const Vec *data)
{
    dh_ = dh;
    map1_ = map1;
    map2_ = map2;
    map3_ = map3;
    data_vec_ = data;
    VecGetArray(*data_vec_, &data_);

    unsigned int ndofs = max(dh_->fe<1>()->n_dofs(), max(dh_->fe<2>()->n_dofs(), dh_->fe<3>()->n_dofs()));
    dof_indices = new unsigned int[ndofs];
}



/**
 * Returns one value in one given point. ResultType can be used to avoid some costly calculation if the result is trivial.
 */
template <int spacedim, class Value>
typename Value::return_type const & FieldFE<spacedim, Value>::value(const Point &p, const ElementAccessor<spacedim> &elm)
{
	Point p_rel = p - elm.element()->node[0]->point();
	DOFHandlerBase::CellIterator cell = dh_->mesh()->element(elm.idx());

	if (elm.dim() == 1) {
		arma::mat::fixed<3,1> m1 = elm.element()->node[1]->point() - elm.element()->node[0]->point();
		arma::mat::fixed<1,3> im1 = pinv(m1);

		Quadrature<1> q1(1);
		q1.set_point(0, im1*p_rel);

		FEValues<1,3> fe_values1(*map1_, q1, *dh_->fe<1>(), update_values);
		fe_values1.reinit(cell);

		dh_->get_dof_indices(cell, dof_indices);

		if (dh_->fe<1>()->is_scalar()) {
			double value = 0;
			for (int i=0; i<dh_->fe<1>()->n_dofs(); i++)
				value += data_[dof_indices[i]]*fe_values1.shape_value(i, 0);
			this->value_(0,0) = value;
		}
		else {
			arma::vec3 value;
			value.zeros();
			for (int i=0; i<dh_->fe<1>()->n_dofs(); i++)
				value += data_[dof_indices[i]]*fe_values1.shape_vector(i, 0);
			for (int i=0; i<3; i++)
				this->value_(i,0) = value(i);
		}
	}
	else if (elm.dim() == 2) {
		arma::mat::fixed<3,2> m2;
		m2.col(0) = elm.element()->node[1]->point() - elm.element()->node[0]->point();
		m2.col(1) = elm.element()->node[2]->point() - elm.element()->node[0]->point();
		arma::mat::fixed<2,3> im2 = pinv(m2);

		Quadrature<2> q2(1);
		q2.set_point(0, im2*p_rel);

		FEValues<2,3> fe_values2(*map2_, q2, *dh_->fe<2>(), update_values);
		fe_values2.reinit(cell);

		dh_->get_dof_indices(cell, dof_indices);

		if (dh_->fe<2>()->is_scalar()) {
			double value = 0;
			for (int i=0; i<dh_->fe<2>()->n_dofs(); i++)
				value += data_[dof_indices[i]]*fe_values2.shape_value(i, 0);
			this->value_(0,0) = value;
		}
		else {
			arma::vec3 value;
			value.zeros();
			for (int i=0; i<dh_->fe<2>()->n_dofs(); i++)
				value += data_[dof_indices[i]]*fe_values2.shape_vector(i, 0);
			for (int i=0; i<3; i++)
				this->value_(i,0) = value(i);
		}	}
	else {
		arma::mat33 m3;
		m3.col(0) = elm.element()->node[1]->point() - elm.element()->node[0]->point();
		m3.col(1) = elm.element()->node[2]->point() - elm.element()->node[0]->point();
		m3.col(2) = elm.element()->node[3]->point() - elm.element()->node[0]->point();
		arma::mat33 im3 = inv(m3);

		Quadrature<3> q3(1);
		q3.set_point(0, im3*p_rel);

		FEValues<3,3> fe_values3(*map3_, q3, *dh_->fe<3>(), update_values);
		fe_values3.reinit(cell);

		dh_->get_dof_indices(cell, dof_indices);

		if (dh_->fe<3>()->is_scalar()) {
			double value = 0;
			for (int i=0; i<dh_->fe<3>()->n_dofs(); i++)
				value += data_[dof_indices[i]]*fe_values3.shape_value(i, 0);
			this->value_(0,0) = value;
		}
		else {
			arma::vec3 value;
			value.zeros();
			for (int i=0; i<dh_->fe<3>()->n_dofs(); i++)
				value += data_[dof_indices[i]]*fe_values3.shape_vector(i, 0);
			for (int i=0; i<3; i++)
				this->value_(i,0) = value(i);
		}
	}

    return this->r_value_;
}



/**
 * Returns std::vector of scalar values in several points at once.
 */
template <int spacedim, class Value>
void FieldFE<spacedim, Value>::value_list (const std::vector< Point >  &point_list, const ElementAccessor<spacedim> &elm,
                   std::vector<typename Value::return_type>  &value_list)
{
    ASSERT_EQUAL( point_list.size(), value_list.size() );

    DOFHandlerBase::CellIterator cell = dh_->mesh()->element( elm.idx() );

	if (elm.dim() == 1) {
		arma::mat::fixed<3,1> m1 = elm.element()->node[1]->point() - elm.element()->node[0]->point();
		arma::mat::fixed<1,3> im1 = pinv(m1);

		dh_->get_dof_indices(cell, dof_indices);

		for (int k=0; k<point_list.size(); k++) {
			Quadrature<1> q1(1);
			Point p_rel = point_list[k] - elm.element()->node[0]->point();
			q1.set_point(0, im1*p_rel);

			FEValues<1,3> fe_values1(*map1_, q1, *dh_->fe<1>(), update_values);
			fe_values1.reinit(cell);

			Value envelope(value_list[k]);

			if (dh_->fe<1>()->is_scalar()) {
				double value = 0;
				for (int i=0; i<dh_->fe<1>()->n_dofs(); i++)
					value += data_[dof_indices[i]]*fe_values1.shape_value(i, 0);
				envelope(0,0) = value;
			}
			else {
				arma::vec3 value;
				value.zeros();
				for (int i=0; i<dh_->fe<1>()->n_dofs(); i++)
					value += data_[dof_indices[i]]*fe_values1.shape_vector(i, 0);
				for (int i=0; i<3; i++)
					envelope(i,0) = value(i);
			}
		}
	}
	else if (elm.dim() == 2) {
		arma::mat::fixed<3,2> m2;
		m2.col(0) = elm.element()->node[1]->point() - elm.element()->node[0]->point();
		m2.col(1) = elm.element()->node[2]->point() - elm.element()->node[0]->point();
		arma::mat::fixed<2,3> im2 = pinv(m2);

		dh_->get_dof_indices(cell, dof_indices);

		for (int k=0; k<point_list.size(); k++) {
			Quadrature<2> q2(1);
			Point p_rel = point_list[k] - elm.element()->node[0]->point();
			q2.set_point(0, im2*p_rel);

			FEValues<2,3> fe_values2(*map2_, q2, *dh_->fe<2>(), update_values);
			fe_values2.reinit(cell);

			Value envelope(value_list[k]);

			if (dh_->fe<2>()->is_scalar()) {
				double value = 0;
				for (int i=0; i<dh_->fe<2>()->n_dofs(); i++)
					value += data_[dof_indices[i]]*fe_values2.shape_value(i, 0);
				envelope(0,0) = value;
			}
			else {
				arma::vec3 value;
				value.zeros();
				for (int i=0; i<dh_->fe<2>()->n_dofs(); i++)
					value += data_[dof_indices[i]]*fe_values2.shape_vector(i, 0);
				for (int i=0; i<3; i++)
					envelope(i,0) = value(i);
			}
		}
	}
	else {
		arma::mat33 m3;
		m3.col(0) = elm.element()->node[1]->point() - elm.element()->node[0]->point();
		m3.col(1) = elm.element()->node[2]->point() - elm.element()->node[0]->point();
		m3.col(2) = elm.element()->node[3]->point() - elm.element()->node[0]->point();
		arma::mat33 im3 = inv(m3);

		dh_->get_dof_indices(cell, dof_indices);

		for (int k=0; k<point_list.size(); k++) {
			Quadrature<3> q3(1);
			Point p_rel = point_list[k] - elm.element()->node[0]->point();
			q3.set_point(0, im3*p_rel);

			FEValues<3,3> fe_values3(*map3_, q3, *dh_->fe<3>(), update_values);
			fe_values3.reinit(cell);

			Value envelope(value_list[k]);

			if (dh_->fe<3>()->is_scalar()) {
				double value = 0;
				for (int i=0; i<dh_->fe<3>()->n_dofs(); i++)
					value += data_[dof_indices[i]]*fe_values3.shape_value(i, 0);
				envelope(0,0) = value;
			}
			else {
				arma::vec3 value;
				value.zeros();
				for (int i=0; i<dh_->fe<3>()->n_dofs(); i++)
					value += data_[dof_indices[i]]*fe_values3.shape_vector(i, 0);
				for (int i=0; i<3; i++)
					envelope(i,0) = value(i);
			}
		}
	}
}



template <int spacedim, class Value>
FieldFE<spacedim, Value>::~FieldFE()
{
	if (dof_indices != nullptr) delete[] dof_indices;
}



#endif /* FIELD_FE_IMPL_HH_ */
