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
 * @file    field_fe.cc
 * @brief   
 */


#include "fields/field_fe.hh"
#include "fields/field_instances.hh"	// for instantiation macros
#include "input/input_type.hh"
#include "quadrature/quadrature.hh"
#include "fem/fe_values.hh"
#include "fem/finite_element.hh"
#include "fem/mapping_p1.hh"
#include "fem/fe_p.hh"
#include "mesh/reader_instances.hh"




/// Implementation.

namespace it = Input::Type;




FLOW123D_FORCE_LINK_IN_CHILD(field_fe)


template <int spacedim, class Value>
const Input::Type::Record & FieldFE<spacedim, Value>::get_input_type()
{
    return it::Record("FieldFE", FieldAlgorithmBase<spacedim,Value>::template_name()+" Field given by finite element approximation.")
        .derive_from(FieldAlgorithmBase<spacedim, Value>::get_input_type())
        .copy_keys(FieldAlgorithmBase<spacedim, Value>::get_field_algo_common_keys())
        .declare_key("mesh_data_file", IT::FileName::input(), IT::Default::obligatory(),
                "GMSH mesh with data. Can be different from actual computational mesh.")
        .declare_key("field_name", IT::String(), IT::Default::obligatory(),
                "The values of the Field are read from the ```$ElementData``` section with field name given by this key.")
        .close();
}

template <int spacedim, class Value>
const int FieldFE<spacedim, Value>::registrar =
		Input::register_class< FieldFE<spacedim, Value>, unsigned int >("FieldFE") +
		FieldFE<spacedim, Value>::get_input_type().size();



template <int spacedim, class Value>
FieldFE<spacedim, Value>::FieldFE( unsigned int n_comp)
: FieldAlgorithmBase<spacedim, Value>(n_comp),
  dh_(nullptr),
  data_vec_(nullptr),
  dof_indices(nullptr),
  map1_(nullptr),
  map2_(nullptr),
  map3_(nullptr),
  field_name_("")
{}


template <int spacedim, class Value>
void FieldFE<spacedim, Value>::set_fe_data(DOFHandlerMultiDim *dh,
		Mapping<1,3> *map1,
		Mapping<2,3> *map2,
		Mapping<3,3> *map3,
		VectorSeqDouble *data)
{
    dh_ = dh;
    map1_ = map1;
    map2_ = map2;
    map3_ = map3;
    data_vec_ = data;

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

		dh_->get_loc_dof_indices(cell, dof_indices);

		if (dh_->fe<1>()->is_scalar()) {
			double value = 0;
			for (unsigned int i=0; i<dh_->fe<1>()->n_dofs(); i++)
				value += (*data_vec_)[dof_indices[i]]*fe_values1.shape_value(i, 0);
			this->value_(0,0) = value;
		}
		else {
			arma::vec3 value;
			value.zeros();
			for (unsigned int i=0; i<dh_->fe<1>()->n_dofs(); i++)
				value += (*data_vec_)[dof_indices[i]]*fe_values1.shape_vector(i, 0);
			for (unsigned int i=0; i<3; i++)
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

		dh_->get_loc_dof_indices(cell, dof_indices);

		if (dh_->fe<2>()->is_scalar()) {
			double value = 0;
			for (unsigned int i=0; i<dh_->fe<2>()->n_dofs(); i++)
				value += (*data_vec_)[dof_indices[i]]*fe_values2.shape_value(i, 0);
			this->value_(0,0) = value;
		}
		else {
			arma::vec3 value;
			value.zeros();
			for (unsigned int i=0; i<dh_->fe<2>()->n_dofs(); i++)
				value += (*data_vec_)[dof_indices[i]]*fe_values2.shape_vector(i, 0);
			for (unsigned int i=0; i<3; i++)
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

		dh_->get_loc_dof_indices(cell, dof_indices);

		if (dh_->fe<3>()->is_scalar()) {
			double value = 0;
			for (unsigned int i=0; i<dh_->fe<3>()->n_dofs(); i++)
				value += (*data_vec_)[dof_indices[i]]*fe_values3.shape_value(i, 0);
			this->value_(0,0) = value;
		}
		else {
			arma::vec3 value;
			value.zeros();
			for (unsigned int i=0; i<dh_->fe<3>()->n_dofs(); i++)
				value += (*data_vec_)[dof_indices[i]]*fe_values3.shape_vector(i, 0);
			for (unsigned int i=0; i<3; i++)
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
	OLD_ASSERT_EQUAL( point_list.size(), value_list.size() );

    DOFHandlerBase::CellIterator cell = dh_->mesh()->element( elm.idx() );

	if (elm.dim() == 1) {
		arma::mat::fixed<3,1> m1 = elm.element()->node[1]->point() - elm.element()->node[0]->point();
		arma::mat::fixed<1,3> im1 = pinv(m1);

		dh_->get_loc_dof_indices(cell, dof_indices);

		for (unsigned int k=0; k<point_list.size(); k++) {
			Quadrature<1> q1(1);
			Point p_rel = point_list[k] - elm.element()->node[0]->point();
			q1.set_point(0, im1*p_rel);

			FEValues<1,3> fe_values1(*map1_, q1, *dh_->fe<1>(), update_values);
			fe_values1.reinit(cell);

			Value envelope(value_list[k]);

			if (dh_->fe<1>()->is_scalar()) {
				double value = 0;
				for (unsigned int i=0; i<dh_->fe<1>()->n_dofs(); i++)
					value += (*data_vec_)[dof_indices[i]]*fe_values1.shape_value(i, 0);
				envelope(0,0) = value;
			}
			else {
				arma::vec3 value;
				value.zeros();
				for (unsigned int i=0; i<dh_->fe<1>()->n_dofs(); i++)
					value += (*data_vec_)[dof_indices[i]]*fe_values1.shape_vector(i, 0);
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

		dh_->get_loc_dof_indices(cell, dof_indices);

		for (unsigned int k=0; k<point_list.size(); k++) {
			Quadrature<2> q2(1);
			Point p_rel = point_list[k] - elm.element()->node[0]->point();
			q2.set_point(0, im2*p_rel);

			FEValues<2,3> fe_values2(*map2_, q2, *dh_->fe<2>(), update_values);
			fe_values2.reinit(cell);

			Value envelope(value_list[k]);

			if (dh_->fe<2>()->is_scalar()) {
				double value = 0;
				for (unsigned int i=0; i<dh_->fe<2>()->n_dofs(); i++)
					value += (*data_vec_)[dof_indices[i]]*fe_values2.shape_value(i, 0);
				envelope(0,0) = value;
			}
			else {
				arma::vec3 value;
				value.zeros();
				for (unsigned int i=0; i<dh_->fe<2>()->n_dofs(); i++)
					value += (*data_vec_)[dof_indices[i]]*fe_values2.shape_vector(i, 0);
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

		dh_->get_loc_dof_indices(cell, dof_indices);

		for (unsigned int k=0; k<point_list.size(); k++) {
			Quadrature<3> q3(1);
			Point p_rel = point_list[k] - elm.element()->node[0]->point();
			q3.set_point(0, im3*p_rel);

			FEValues<3,3> fe_values3(*map3_, q3, *dh_->fe<3>(), update_values);
			fe_values3.reinit(cell);

			Value envelope(value_list[k]);

			if (dh_->fe<3>()->is_scalar()) {
				double value = 0;
				for (unsigned int i=0; i<dh_->fe<3>()->n_dofs(); i++)
					value += (*data_vec_)[dof_indices[i]]*fe_values3.shape_value(i, 0);
				envelope(0,0) = value;
			}
			else {
				arma::vec3 value;
				value.zeros();
				for (unsigned int i=0; i<dh_->fe<3>()->n_dofs(); i++)
					value += (*data_vec_)[dof_indices[i]]*fe_values3.shape_vector(i, 0);
				for (int i=0; i<3; i++)
					envelope(i,0) = value(i);
			}
		}
	}
}



template <int spacedim, class Value>
void FieldFE<spacedim, Value>::init_from_input(const Input::Record &rec, const struct FieldAlgoBaseInitData& init_data) {
	this->init_unit_conversion_coefficient(rec, init_data);
	flags_ = init_data.flags_;


	// read mesh, create tree
    {
       reader_file_ = FilePath( rec.val<FilePath>("mesh_data_file") );
       source_mesh_ = ReaderInstance::get_mesh(reader_file_);
       source_mesh_->get_bih_tree(); // only create BIH tree
	   // no call to mesh->setup_topology, we need only elements, no connectivity
    }

	field_name_ = rec.val<std::string>("field_name");
}



template <int spacedim, class Value>
void FieldFE<spacedim, Value>::set_mesh(const Mesh *mesh, bool boundary_domain) {
	// Mesh can be set only for field initialized from input.
	if ( flags_.match(FieldFlag::equation_input) && flags_.match(FieldFlag::declare_input) ) {
		ASSERT(field_name_ != "").error("Uninitialized FieldFE, did you call init_from_input()?\n");
		this->make_dof_handler(mesh);
	}
}



template <int spacedim, class Value>
void FieldFE<spacedim, Value>::make_dof_handler(const Mesh *mesh) {
	// temporary solution - these objects will be set through FieldCommon
	map1_ = new MappingP1<1,3>();
	map2_ = new MappingP1<2,3>();
	map3_ = new MappingP1<3,3>();
	fe1_ = new FE_P_disc<0,1,3>();
	fe2_ = new FE_P_disc<0,2,3>();
	fe3_ = new FE_P_disc<0,3,3>();

	dh_ = new DOFHandlerMultiDim( const_cast<Mesh &>(*mesh) );
	dh_->distribute_dofs(*fe1_, *fe2_, *fe3_);
    unsigned int ndofs = max(dh_->fe<1>()->n_dofs(), max(dh_->fe<2>()->n_dofs(), dh_->fe<3>()->n_dofs()));
    dof_indices = new unsigned int[ndofs];

    // allocate data_vec_
	unsigned int data_size = dh_->n_global_dofs();
	data_vec_ = new VectorSeqDouble();
	data_vec_->resize(data_size);
}



template <int spacedim, class Value>
bool FieldFE<spacedim, Value>::set_time(const TimeStep &time) {
	// Time can be set only for field initialized from input.
	if ( flags_.match(FieldFlag::equation_input) && flags_.match(FieldFlag::declare_input) ) {
	    ASSERT(field_name_ != "").error("Uninitialized FieldFE, did you call init_from_input()?\n");
		ASSERT_PTR(dh_)(field_name_).error("Null target mesh pointer of finite element field, did you call set_mesh()?\n");
		if ( reader_file_ == FilePath() ) return false;

		GMSH_DataHeader search_header;
		search_header.actual = false;
		search_header.field_name = field_name_;
		search_header.n_components = this->value_.n_rows() * this->value_.n_cols();
		search_header.n_entities = source_mesh_->element.size();
		search_header.time = time.end();

		bool boundary_domain_ = false;
		std::vector<double> data_vec = *(ReaderInstance::get_reader(reader_file_)->template get_element_data<double>(search_header,
				source_mesh_->elements_id_maps(boundary_domain_), this->component_idx_));
		std::vector<double> sum_val(4);
		std::vector<unsigned int> elem_count(4);

		FOR_ELEMENTS( dh_->mesh(), ele ) {
			searched_elements_.clear();
			source_mesh_->get_bih_tree().find_point(ele->centre(), searched_elements_);
			std::fill(sum_val.begin(), sum_val.end(), 0.0);
			std::fill(elem_count.begin(), elem_count.end(), 0);
			for (std::vector<unsigned int>::iterator it = searched_elements_.begin(); it!=searched_elements_.end(); it++) {
				ElementFullIter elm = source_mesh_->element( *it );
				arma::mat map = elm->element_map();
				arma::vec projection = elm->project_point(ele->centre(), map);
				if (projection.min() >= -BoundingBox::epsilon) {
					// projection in element
					arma::vec3 projection_point = map.col(elm->dim()) + map.cols(0, elm->dim()-1) * projection.rows(0, elm->dim()-1);
					if ( arma::norm(projection_point - ele->centre(), "inf") < 2*BoundingBox::epsilon ) {
						// point on the element
						sum_val[elm->dim()] += data_vec[*it];
						++elem_count[elm->dim()];
					}
				}
			}
			unsigned int dim = ele->dim();
			double elem_value = 0.0;
			do {
				if (elem_count[dim] > 0) {
					elem_value = sum_val[dim] / elem_count[dim];
					break;
				}
				++dim;
			} while (dim<4);

			dh_->get_loc_dof_indices( ele, dof_indices);
			ASSERT_LT_DBG( dof_indices[0], data_vec_->size());
			(*data_vec_)[dof_indices[0]] = elem_value * this->unit_conversion_coefficient_;
		}


		return search_header.actual;
	} else return false;

}



template <int spacedim, class Value>
FieldFE<spacedim, Value>::~FieldFE()
{
	if (dof_indices != nullptr) delete[] dof_indices;
}


// Instantiations of FieldFE
INSTANCE_ALL(FieldFE)
