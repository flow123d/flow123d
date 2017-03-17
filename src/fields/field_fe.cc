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
  data_vec_(nullptr),
  dof_indices(nullptr),
  field_name_("")
{}


template <int spacedim, class Value>
void FieldFE<spacedim, Value>::set_fe_data(std::shared_ptr<DOFHandlerMultiDim> dh,
		Mapping<1,3> *map1,
		Mapping<2,3> *map2,
		Mapping<3,3> *map3,
		VectorSeqDouble *data)
{
    dh_ = dh;
    data_vec_ = data;

    unsigned int ndofs = dh_->max_elem_dofs();
    dof_indices = new unsigned int[ndofs];

    // initialization data of value handlers
	FEValueInitData init_data;
	init_data.dh = dh_;
	init_data.data_vec = data_vec_;
	init_data.ndofs = ndofs;
	init_data.n_comp = this->n_comp();

	// initialize value handler objects
	value_handler1_.initialize(init_data, map1);
	value_handler2_.initialize(init_data, map2);
	value_handler3_.initialize(init_data, map3);
}



/**
 * Returns one value in one given point. ResultType can be used to avoid some costly calculation if the result is trivial.
 */
template <int spacedim, class Value>
typename Value::return_type const & FieldFE<spacedim, Value>::value(const Point &p, const ElementAccessor<spacedim> &elm)
{
	switch (elm.dim()) {
	case 1:
		return value_handler1_.value(p, elm);
	case 2:
		return value_handler2_.value(p, elm);
	case 3:
		return value_handler3_.value(p, elm);
	default:
		ASSERT(false).error("Invalid element dimension!");
	}

    return this->r_value_;
}



/**
 * Returns std::vector of scalar values in several points at once.
 */
template <int spacedim, class Value>
void FieldFE<spacedim, Value>::value_list (const std::vector< Point > &point_list, const ElementAccessor<spacedim> &elm,
                   std::vector<typename Value::return_type> &value_list)
{
	ASSERT_EQ( point_list.size(), value_list.size() ).error();

	switch (elm.dim()) {
	case 1:
		value_handler1_.value_list(point_list, elm, value_list);
		break;
	case 2:
		value_handler2_.value_list(point_list, elm, value_list);
		break;
	case 3:
		value_handler3_.value_list(point_list, elm, value_list);
		break;
	default:
		ASSERT(false).error("Invalid element dimension!");
	}
}



template <int spacedim, class Value>
void FieldFE<spacedim, Value>::init_from_input(const Input::Record &rec, const struct FieldAlgoBaseInitData& init_data) {
	this->init_unit_conversion_coefficient(rec, init_data);
	flags_ = init_data.flags_;


	// read data from input record
    reader_file_ = FilePath( rec.val<FilePath>("mesh_data_file") );
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
	fe1_ = new FE_P_disc<0,1,3>();
	fe2_ = new FE_P_disc<0,2,3>();
	fe3_ = new FE_P_disc<0,3,3>();

	dh_ = std::make_shared<DOFHandlerMultiDim>( const_cast<Mesh &>(*mesh) );
    std::shared_ptr<DiscreteSpace> ds = std::make_shared<EqualOrderDiscreteSpace>( &const_cast<Mesh &>(*mesh), fe1_, fe2_, fe3_);
	dh_->distribute_dofs(ds);
    unsigned int ndofs = dh_->max_elem_dofs();
    dof_indices = new unsigned int[ndofs];

    // allocate data_vec_
	unsigned int data_size = dh_->n_global_dofs();
	data_vec_ = new VectorSeqDouble();
	data_vec_->resize(data_size);

	// initialization data of value handlers
	FEValueInitData init_data;
	init_data.dh = dh_;
	init_data.data_vec = data_vec_;
	init_data.ndofs = ndofs;
	init_data.n_comp = this->n_comp();

	// initialize value handler objects
	value_handler1_.initialize(init_data);
	value_handler2_.initialize(init_data);
	value_handler3_.initialize(init_data);
}



template <int spacedim, class Value>
bool FieldFE<spacedim, Value>::set_time(const TimeStep &time) {
	// Time can be set only for field initialized from input.
	if ( flags_.match(FieldFlag::equation_input) && flags_.match(FieldFlag::declare_input) ) {
	    ASSERT(field_name_ != "").error("Uninitialized FieldFE, did you call init_from_input()?\n");
		ASSERT_PTR(dh_)(field_name_).error("Null target mesh pointer of finite element field, did you call set_mesh()?\n");
		if ( reader_file_ == FilePath() ) return false;

		std::shared_ptr<Mesh> source_mesh = ReaderInstance::get_mesh(reader_file_);

		GMSH_DataHeader search_header;
		search_header.actual = false;
		search_header.field_name = field_name_;
		search_header.n_components = this->value_.n_rows() * this->value_.n_cols();
		search_header.n_entities = source_mesh->element.size();
		search_header.time = time.end();

		bool boundary_domain_ = false;
		auto data_vec = ReaderInstance::get_reader(reader_file_)->template get_element_data<double>(search_header,
				source_mesh->elements_id_maps(boundary_domain_), this->component_idx_);
		this->interpolate(data_vec);

		return search_header.actual;
	} else return false;

}


template <int spacedim, class Value>
void FieldFE<spacedim, Value>::interpolate(ElementDataCache<double>::ComponentDataPtr data_vec)
{
	std::shared_ptr<Mesh> source_mesh = ReaderInstance::get_mesh(reader_file_);
	std::vector<double> sum_val(4);
	std::vector<unsigned int> elem_count(4);
	std::vector<unsigned int> searched_elements; // stored suspect elements in calculating the intersection

	FOR_ELEMENTS( dh_->mesh(), ele ) {
		searched_elements.clear();
		source_mesh->get_bih_tree().find_point(ele->centre(), searched_elements);
		std::fill(sum_val.begin(), sum_val.end(), 0.0);
		std::fill(elem_count.begin(), elem_count.end(), 0);
		for (std::vector<unsigned int>::iterator it = searched_elements.begin(); it!=searched_elements.end(); it++) {
			ElementFullIter elm = source_mesh->element( *it );
			arma::mat map = elm->element_map();
			arma::vec projection = elm->project_point(ele->centre(), map);
			if (projection.min() >= -BoundingBox::epsilon) {
				// projection in element
				arma::vec3 projection_point = map.col(elm->dim()) + map.cols(0, elm->dim()-1) * projection.rows(0, elm->dim()-1);
				if ( arma::norm(projection_point - ele->centre(), "inf") < 2*BoundingBox::epsilon ) {
					// point on the element
					sum_val[elm->dim()] += (*data_vec)[*it];
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
}



template <int spacedim, class Value>
FieldFE<spacedim, Value>::~FieldFE()
{
	if (dof_indices != nullptr) delete[] dof_indices;
}


// Instantiations of FieldFE
INSTANCE_ALL(FieldFE)
