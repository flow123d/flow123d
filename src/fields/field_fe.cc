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


#include <limits>

#include "fields/field_fe.hh"
#include "fields/vec_seq_double.hh"
#include "fields/field_instances.hh"	// for instantiation macros
#include "fields/fe_value_handler.hh"
#include "input/input_type.hh"
#include "fem/fe_p.hh"
#include "fem/fe_system.hh"
#include "fem/dh_cell_accessor.hh"
#include "io/reader_cache.hh"
#include "io/msh_gmshreader.h"
#include "mesh/accessors.hh"
#include "mesh/range_wrapper.hh"




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
        .declare_key("input_discretization", FieldFE<spacedim, Value>::get_disc_selection_input_type(), IT::Default::optional(),
                "Section where to find the field, some sections are specific to file format \n"
        		"point_data/node_data, cell_data/element_data, -/element_node_data, native/-.\n"
        		"If not given by user we try to find the field in all sections, but report error \n"
        		"if it is found in more the one section.")
        .declare_key("field_name", IT::String(), IT::Default::obligatory(),
                "The values of the Field are read from the ```$ElementData``` section with field name given by this key.")
        .declare_key("default_value", IT::Double(), IT::Default::optional(),
                "Allow set default value of elements that have not listed values in mesh data file.")
        .declare_key("time_unit", IT::String(), IT::Default::read_time("Common unit of TimeGovernor."),
                "Definition of unit of all times defined in mesh data file.")
        .declare_key("read_time_shift", TimeGovernor::get_input_time_type(), IT::Default("0.0"),
                "Allow set time shift of field data read from the mesh data file. For time 't', field descriptor with time 'T', "
                "time shift 'S' and if 't > T', we read time frame 't + S'.")
        .close();
}

template <int spacedim, class Value>
const Input::Type::Selection & FieldFE<spacedim, Value>::get_disc_selection_input_type()
{
	return it::Selection("FE_discretization",
			"Specify the section in mesh input file where field data is listed.\nSome sections are specific to file format.")
		.add_value(OutputTime::DiscreteSpace::NODE_DATA, "node_data", "point_data (VTK) / node_data (GMSH)")
		.add_value(OutputTime::DiscreteSpace::ELEM_DATA, "element_data", "cell_data (VTK) / element_data (GMSH)")
		.add_value(OutputTime::DiscreteSpace::CORNER_DATA, "element_node_data", "element_node_data (only for GMSH)")
		.add_value(OutputTime::DiscreteSpace::NATIVE_DATA, "native_data", "native_data (only for VTK)")
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
  field_name_(""),
  has_compatible_mesh_(false)
{
	this->is_constant_in_space_ = false;
}


template <int spacedim, class Value>
void FieldFE<spacedim, Value>::set_fe_data(std::shared_ptr<DOFHandlerMultiDim> dh,
		MappingP1<1,3> *map1,
		MappingP1<2,3> *map2,
		MappingP1<3,3> *map3,
		VectorSeqDouble *data)
{
    dh_ = dh;
    data_vec_ = data;
    reinit_fe_data(map1, map2, map3);
}



template <int spacedim, class Value>
void FieldFE<spacedim, Value>::set_native_dh(std::shared_ptr<DOFHandlerMultiDim> dh)
{
    dh_ = dh;
    reinit_fe_data(nullptr, nullptr, nullptr);
}



template <int spacedim, class Value>
void FieldFE<spacedim, Value>::reinit_fe_data(MappingP1<1,3> *map1,
		MappingP1<2,3> *map2,
		MappingP1<3,3> *map3)
{
	//ASSERT_EQ(dh_->n_global_dofs(), data_vec_->size());

    unsigned int ndofs = dh_->max_elem_dofs();
    dof_indices_.resize(ndofs);

    // initialization data of value handlers
	FEValueInitData init_data;
	init_data.dh = dh_;
	init_data.data_vec = data_vec_;
	init_data.ndofs = ndofs;
	init_data.n_comp = this->n_comp();

	// initialize value handler objects
	value_handler0_.initialize(init_data);
	value_handler1_.initialize(init_data, map1);
	value_handler2_.initialize(init_data, map2);
	value_handler3_.initialize(init_data, map3);

	// set discretization
	discretization_ = OutputTime::DiscreteSpace::UNDEFINED;
}



/**
 * Returns one value in one given point. ResultType can be used to avoid some costly calculation if the result is trivial.
 */
template <int spacedim, class Value>
typename Value::return_type const & FieldFE<spacedim, Value>::value(const Point &p, const ElementAccessor<spacedim> &elm)
{
	switch (elm.dim()) {
	case 0:
		return value_handler0_.value(p, elm);
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
	case 0:
		value_handler0_.value_list(point_list, elm, value_list);
		break;
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
	this->in_rec_ = rec;
	flags_ = init_data.flags_;


	// read data from input record
    reader_file_ = FilePath( rec.val<FilePath>("mesh_data_file") );
	field_name_ = rec.val<std::string>("field_name");
	if (! rec.opt_val<OutputTime::DiscreteSpace>("input_discretization", discretization_) ) {
		discretization_ = OutputTime::DiscreteSpace::UNDEFINED;
	}
    if (! rec.opt_val("default_value", default_value_) ) {
    	default_value_ = numeric_limits<double>::signaling_NaN();
    }
}



template <int spacedim, class Value>
void FieldFE<spacedim, Value>::set_mesh(const Mesh *mesh, bool boundary_domain) {
	// Mesh can be set only for field initialized from input.
	if ( flags_.match(FieldFlag::equation_input) && flags_.match(FieldFlag::declare_input) ) {
		ASSERT(field_name_ != "").error("Uninitialized FieldFE, did you call init_from_input()?\n");
		this->boundary_domain_ = boundary_domain;
		this->make_dof_handler(mesh);
		this->has_compatible_mesh_ = ReaderCache::check_compatible_mesh(reader_file_, const_cast<Mesh &>(*mesh) );
	}
}



template <int spacedim, class Value>
void FieldFE<spacedim, Value>::make_dof_handler(const Mesh *mesh) {
	// temporary solution - these objects will be set through FieldCommon
	switch (this->value_.n_rows() * this->value_.n_cols()) { // by number of components
		case 1: { // scalar
			fe0_ = new FE_P_disc<0>(0);
			fe1_ = new FE_P_disc<1>(0);
			fe2_ = new FE_P_disc<2>(0);
			fe3_ = new FE_P_disc<3>(0);
			break;
		}
		case 3: { // vector
			std::shared_ptr< FiniteElement<0> > fe0_ptr = std::make_shared< FE_P_disc<0> >(0);
			std::shared_ptr< FiniteElement<1> > fe1_ptr = std::make_shared< FE_P_disc<1> >(0);
			std::shared_ptr< FiniteElement<2> > fe2_ptr = std::make_shared< FE_P_disc<2> >(0);
			std::shared_ptr< FiniteElement<3> > fe3_ptr = std::make_shared< FE_P_disc<3> >(0);
			fe0_ = new FESystem<0>(fe0_ptr, FEType::FEVector, 3);
			fe1_ = new FESystem<1>(fe1_ptr, FEType::FEVector, 3);
			fe2_ = new FESystem<2>(fe2_ptr, FEType::FEVector, 3);
			fe3_ = new FESystem<3>(fe3_ptr, FEType::FEVector, 3);
			break;
		}
		case 9: { // tensor
			std::shared_ptr< FiniteElement<0> > fe0_ptr = std::make_shared< FE_P_disc<0> >(0);
			std::shared_ptr< FiniteElement<1> > fe1_ptr = std::make_shared< FE_P_disc<1> >(0);
			std::shared_ptr< FiniteElement<2> > fe2_ptr = std::make_shared< FE_P_disc<2> >(0);
			std::shared_ptr< FiniteElement<3> > fe3_ptr = std::make_shared< FE_P_disc<3> >(0);
			fe0_ = new FESystem<0>(fe0_ptr, FEType::FETensor, 9);
			fe1_ = new FESystem<1>(fe1_ptr, FEType::FETensor, 9);
			fe2_ = new FESystem<2>(fe2_ptr, FEType::FETensor, 9);
			fe3_ = new FESystem<3>(fe3_ptr, FEType::FETensor, 9);
			break;
		}
		default:
			ASSERT(false).error("Should not happen!\n");
	}

	dh_ = std::make_shared<DOFHandlerMultiDim>( const_cast<Mesh &>(*mesh) );
    std::shared_ptr<DiscreteSpace> ds = std::make_shared<EqualOrderDiscreteSpace>( &const_cast<Mesh &>(*mesh), fe0_, fe1_, fe2_, fe3_);
	dh_->distribute_dofs(ds, true);
    unsigned int ndofs = dh_->max_elem_dofs();
    dof_indices_.resize(ndofs);

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
	value_handler0_.initialize(init_data);
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

		unsigned int n_components = this->value_.n_rows() * this->value_.n_cols();
		double time_unit_coef = time.read_coef(in_rec_.find<string>("time_unit"));
		double time_shift = time.read_time( in_rec_.find<Input::Tuple>("read_time_shift") );
		double read_time = (time.end()+time_shift) / time_unit_coef;
		BaseMeshReader::HeaderQuery header_query(field_name_, read_time, this->discretization_, dh_->hash());
		ReaderCache::get_reader(reader_file_)->find_header(header_query);
		// TODO: use default and check NaN values in data_vec

		unsigned int n_entities;
		bool is_native = (header_query.discretization == OutputTime::DiscreteSpace::NATIVE_DATA);
		if (is_native || this->has_compatible_mesh_) {
			n_entities = dh_->mesh()->n_elements();
		} else {
			n_entities = ReaderCache::get_mesh(reader_file_)->n_elements();
		}
		auto data_vec = ReaderCache::get_reader(reader_file_)->template get_element_data<double>(n_entities, n_components,
				this->boundary_domain_, this->component_idx_);
		CheckResult checked_data = ReaderCache::get_reader(reader_file_)->scale_and_check_limits(field_name_,
				this->unit_conversion_coefficient_, default_value_);


	    if (checked_data == CheckResult::not_a_number) {
	        THROW( ExcUndefElementValue() << EI_Field(field_name_) );
	    }

		if (is_native || this->has_compatible_mesh_) {
			this->calculate_native_values(data_vec);
		} else {
			this->interpolate(data_vec);
		}

		return true;
	} else return false;

}


template <int spacedim, class Value>
void FieldFE<spacedim, Value>::interpolate(ElementDataCache<double>::ComponentDataPtr data_vec)
{
	ASSERT(!this->boundary_domain_)(field_name_).error("Interpolation of boundary FieldFE is not supported yet.\n");
	std::shared_ptr<Mesh> source_mesh = ReaderCache::get_mesh(reader_file_);
	std::vector<double> sum_val(4);
	std::vector<unsigned int> elem_count(4);
	std::vector<unsigned int> searched_elements; // stored suspect elements in calculating the intersection

	for (auto ele : dh_->mesh()->elements_range()) {
		searched_elements.clear();
		source_mesh->get_bih_tree().find_point(ele.centre(), searched_elements);
		std::fill(sum_val.begin(), sum_val.end(), 0.0);
		std::fill(elem_count.begin(), elem_count.end(), 0);
		for (std::vector<unsigned int>::iterator it = searched_elements.begin(); it!=searched_elements.end(); it++) {
			ElementAccessor<3> elm = source_mesh->element_accessor(*it);
			bool contains=false;
			switch (elm->dim()) {
			case 0:
				contains = arma::norm(elm.node(0)->point()-ele.centre(), 2) < 4*std::numeric_limits<double>::epsilon();
				break;
			case 1:
				contains = value_handler1_.get_mapping()->contains_point(ele.centre(), elm);
				break;
			case 2:
				contains = value_handler2_.get_mapping()->contains_point(ele.centre(), elm);
				break;
			case 3:
				contains = value_handler3_.get_mapping()->contains_point(ele.centre(), elm);
				break;
			default:
				ASSERT(false).error("Invalid element dimension!");
			}
			if (contains) {
				// projection point in element
				sum_val[elm->dim()] += (*data_vec)[*it];
				++elem_count[elm->dim()];
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

		dh_->get_dof_indices( ele, dof_indices_);
		ASSERT_LT_DBG( dof_indices_[0], (int)data_vec_->size());
		(*data_vec_)[dof_indices_[0]] = elem_value * this->unit_conversion_coefficient_;
	}
}


template <int spacedim, class Value>
void FieldFE<spacedim, Value>::calculate_native_values(ElementDataCache<double>::ComponentDataPtr data_cache)
{
	// Same algorithm as in output of Node_data. Possibly code reuse.
	unsigned int dof_size, data_vec_i;
	std::vector<unsigned int> count_vector(data_vec_->size(), 0);
	data_vec_->fill(0.0);
	VectorSeqDouble::VectorSeq data_vector = data_vec_->get_data_ptr();

	int rank;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);

	if (rank == 0) {
		// iterate through elements, assembly global vector and count number of writes
		for (auto ele : dh_->mesh()->elements_range()) { // remove special case for rank == 0 - necessary for correct output
			dof_size = dh_->get_dof_indices( ele, dof_indices_ );
			data_vec_i = ele.idx() * dof_indices_.size();
			for (unsigned int i=0; i<dof_size; ++i, ++data_vec_i) {
				(*data_vector)[ dof_indices_[i] ] += (*data_cache)[data_vec_i];
				++count_vector[ dof_indices_[i] ];
			}
		}
	} else {
		// iterate through cells, assembly global vector and count number of writes
		for (auto cell : dh_->own_range()) {
			dof_size = cell.get_dof_indices(dof_indices_);
			data_vec_i = cell.element_idx() * dof_indices_.size();
			for (unsigned int i=0; i<dof_size; ++i, ++data_vec_i) {
				(*data_vector)[ dof_indices_[i] ] += (*data_cache)[data_vec_i];
				++count_vector[ dof_indices_[i] ];
			}
		}
	}

	// compute averages of values
	for (unsigned int i=0; i<data_vec_->size(); ++i) {
		if (count_vector[i]>0) (*data_vector)[i] /= count_vector[i];
	}
}


template <int spacedim, class Value>
void FieldFE<spacedim, Value>::fill_data_to_cache(ElementDataCache<double> &output_data_cache) {
	ASSERT_EQ(output_data_cache.n_values() * output_data_cache.n_elem(), dh_->n_global_dofs()).error();
	ASSERT_EQ(output_data_cache.n_elem(), dof_indices_.size()).error();
	double loc_values[output_data_cache.n_elem()];
	unsigned int i, dof_filled_size;

	VectorSeqDouble::VectorSeq data_vec = data_vec_->get_data_ptr();
	for (auto ele : dh_->mesh()->elements_range()) {
		dof_filled_size = dh_->get_dof_indices( ele, dof_indices_);
		for (i=0; i<dof_filled_size; ++i) loc_values[i] = (*data_vec)[ dof_indices_[0] ];
		for ( ; i<output_data_cache.n_elem(); ++i) loc_values[i] = numeric_limits<double>::signaling_NaN();
		output_data_cache.store_value( ele.idx(), loc_values );
	}

	output_data_cache.set_dof_handler_hash( dh_->hash() );
}



template <int spacedim, class Value>
inline unsigned int FieldFE<spacedim, Value>::data_size() const {
	return data_vec_->size();
}



template <int spacedim, class Value>
FieldFE<spacedim, Value>::~FieldFE()
{}


// Instantiations of FieldFE
INSTANCE_ALL(FieldFE)
