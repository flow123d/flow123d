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
#include "la/vector_mpi.hh"
#include "fields/field_instances.hh"	// for instantiation macros
#include "fields/fe_value_handler.hh"
#include "input/input_type.hh"
#include "fem/fe_p.hh"
#include "fem/fe_system.hh"
#include "fem/dh_cell_accessor.hh"
#include "fem/mapping_p1.hh"
#include "io/reader_cache.hh"
#include "io/msh_gmshreader.h"
#include "mesh/accessors.hh"
#include "mesh/range_wrapper.hh"
#include "mesh/bc_mesh.hh"
#include "quadrature/quadrature_lib.hh"

#include "system/sys_profiler.hh"
#include "intersection/intersection_aux.hh"
#include "intersection/intersection_local.hh"
#include "intersection/compute_intersection.hh"




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
                "Section where to find the field.\n Some sections are specific to file format: "
        		"point_data/node_data, cell_data/element_data, -/element_node_data, native/-.\n"
        		"If not given by a user, we try to find the field in all sections, but we report an error "
        		"if it is found in more than one section.")
        .declare_key("field_name", IT::String(), IT::Default::obligatory(),
                "The values of the Field are read from the ```$ElementData``` section with field name given by this key.")
        .declare_key("default_value", IT::Double(), IT::Default::optional(),
                "Default value is set on elements which values have not been listed in the mesh data file.")
        .declare_key("time_unit", IT::String(), IT::Default::read_time("Common unit of TimeGovernor."),
                "Definition of the unit of all times defined in the mesh data file.")
        .declare_key("read_time_shift", TimeGovernor::get_input_time_type(), IT::Default("0.0"),
                "This key allows reading field data from the mesh data file shifted in time. Considering the time 't', field descriptor with time 'T', "
                "time shift 'S', then if 't > T', we read the time frame 't + S'.")
        .declare_key("interpolation", FieldFE<spacedim, Value>::get_interp_selection_input_type(),
        		IT::Default("\"equivalent_mesh\""), "Type of interpolation applied to the input spatial data.\n"
        		"The default value 'equivalent_mesh' assumes the data being constant on elements living on the same mesh "
        		"as the computational mesh, but possibly with different numbering. In the case of the same numbering, "
        		"the user can set 'identical_mesh' to omit algorithm for guessing node and element renumbering. "
        		"Alternatively, in case of different input mesh, several interpolation algorithms are available.")
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
const Input::Type::Selection & FieldFE<spacedim, Value>::get_interp_selection_input_type()
{
	return it::Selection("interpolation", "Specify interpolation of the input data from its input mesh to the computational mesh.")
		.add_value(DataInterpolation::identic_msh, "identic_mesh", "Topology and indices of nodes and elements of"
				"the input mesh and the computational mesh are identical. "
				"This interpolation is typically used for GMSH input files containing only the field values without "
				"explicit mesh specification.")
		.add_value(DataInterpolation::equivalent_msh, "equivalent_mesh", "Topologies of the input mesh and the computational mesh "
				"are the same, the node and element numbering may differ. "
				"This interpolation can be used also for VTK input data.") // default value
		.add_value(DataInterpolation::gauss_p0, "P0_gauss", "Topologies of the input mesh and the computational mesh may differ. "
				"Constant values on the elements of the computational mesh are evaluated using the Gaussian quadrature of the fixed order 4, "
				"where the quadrature points and their values are found in the input mesh and input data using the BIH tree search."
				)
		.add_value(DataInterpolation::interp_p0, "P0_intersection", "Topologies of the input mesh and the computational mesh may differ. "
				"Can be applied only for boundary fields. For every (boundary) element of the computational mesh the "
				"intersection with the input mesh is computed. Constant values on the elements of the computational mesh "
				"are evaluated as the weighted average of the (constant) values on the intersecting elements of the input mesh.")
		.close();
}

template <int spacedim, class Value>
const int FieldFE<spacedim, Value>::registrar =
		Input::register_class< FieldFE<spacedim, Value>, unsigned int >("FieldFE") +
		FieldFE<spacedim, Value>::get_input_type().size();



template <int spacedim, class Value>
FieldFE<spacedim, Value>::FieldFE( unsigned int n_comp)
: FieldAlgorithmBase<spacedim, Value>(n_comp),
  field_name_("")
{
	this->is_constant_in_space_ = false;
}


template <int spacedim, class Value>
VectorMPI FieldFE<spacedim, Value>::set_fe_data(std::shared_ptr<DOFHandlerMultiDim> dh,
		unsigned int component_index, VectorMPI dof_values)
{
    dh_ = dh;
    if (dof_values.size()==0) { //create data vector according to dof handler - Warning not tested yet
        data_vec_ = dh_->create_vector();
        data_vec_.zero_entries();
    } else {
        data_vec_ = dof_values;
    }

    unsigned int ndofs = dh_->max_elem_dofs();

    // initialization data of value handlers
	FEValueInitData init_data;
	init_data.dh = dh_;
	init_data.data_vec = data_vec_;
	init_data.ndofs = ndofs;
	init_data.n_comp = this->n_comp();
	init_data.comp_index = component_index;

	// initialize value handler objects
	value_handler0_.initialize(init_data);
	value_handler1_.initialize(init_data);
	value_handler2_.initialize(init_data);
	value_handler3_.initialize(init_data);

	// set discretization
	discretization_ = OutputTime::DiscreteSpace::UNDEFINED;
	interpolation_ = DataInterpolation::equivalent_msh;

	return data_vec_;
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
void FieldFE<spacedim, Value>::value_list (const Armor::array &point_list, const ElementAccessor<spacedim> &elm,
                   std::vector<typename Value::return_type> &value_list)
{
	ASSERT_EQ( point_list.size(), value_list.size() ).error();
	ASSERT_DBG( point_list.n_rows() == spacedim && point_list.n_cols() == 1).error("Invalid point size.\n");

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
	if (! rec.opt_val<DataInterpolation>("interpolation", interpolation_) ) {
		interpolation_ = DataInterpolation::equivalent_msh;
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
		switch (this->interpolation_) {
			case DataInterpolation::identic_msh:
				ReaderCache::get_element_ids(reader_file_, *mesh);
				break;
			case DataInterpolation::equivalent_msh:
				if (!ReaderCache::check_compatible_mesh( reader_file_, const_cast<Mesh &>(*mesh) )) {
					this->interpolation_ = DataInterpolation::gauss_p0;
					WarningOut().fmt("Source mesh of FieldFE '{}' is not compatible with target mesh.\nInterpolation of input data will be changed to 'P0_gauss'.\n",
							field_name_);
				}
				break;
			case DataInterpolation::gauss_p0:
			{
				auto source_mesh = ReaderCache::get_mesh(reader_file_);
				ReaderCache::get_element_ids(reader_file_, *(source_mesh.get()) );
				break;
			}
			case DataInterpolation::interp_p0:
			{
				auto source_msh = ReaderCache::get_mesh(reader_file_);
				ReaderCache::get_element_ids(reader_file_, *(source_msh.get()) );
				if (!boundary_domain) {
					this->interpolation_ = DataInterpolation::gauss_p0;
					WarningOut().fmt("Interpolation 'P0_intersection' of FieldFE '{}' can't be used on bulk region.\nIt will be changed to 'P0_gauss'.\n", field_name_);
				}
				break;
			}
		}
		this->make_dof_handler(mesh);
	}
}



template <int spacedim, class Value>
void FieldFE<spacedim, Value>::fill_boundary_dofs() {
	ASSERT(this->boundary_domain_);

	auto bc_mesh = dh_->mesh()->get_bc_mesh();
	unsigned int n_comp = this->value_.n_rows() * this->value_.n_cols();
	boundary_dofs_ = std::make_shared< std::vector<Idx> >( n_comp * bc_mesh->n_elements() );
	std::vector<Idx> &in_vec = *( boundary_dofs_.get() );
	unsigned int j = 0; // actual index to boundary_dofs_ vector

	for (auto ele : bc_mesh->elements_range()) {
		Idx elm_shift = n_comp * ele.idx();
		for (unsigned int i=0; i<n_comp; ++i, ++j) {
			in_vec[j] = elm_shift + i;
		}
	}

	value_handler0_.set_boundary_dofs_vector(boundary_dofs_);
	value_handler1_.set_boundary_dofs_vector(boundary_dofs_);
	value_handler2_.set_boundary_dofs_vector(boundary_dofs_);
	value_handler3_.set_boundary_dofs_vector(boundary_dofs_);

	data_vec_.resize(boundary_dofs_->size());
}



template <int spacedim, class Value>
void FieldFE<spacedim, Value>::make_dof_handler(const Mesh *mesh) {
	// temporary solution - these objects will be set through FieldCommon
	switch (this->value_.n_rows() * this->value_.n_cols()) { // by number of components
		case 1: { // scalar
			fe_ = MixedPtr<FE_P_disc>(0);
			break;
		}
		case 3: { // vector
			 MixedPtr<FE_P_disc>   fe_base(0) ;
			fe_ = mixed_fe_system(fe_base, FEType::FEVector, 3);
			break;
		}
		case 9: { // tensor
		    MixedPtr<FE_P_disc>   fe_base(0) ;
            fe_ = mixed_fe_system(fe_base, FEType::FETensor, 9);
			break;
		}
		default:
			ASSERT(false).error("Should not happen!\n");
	}

	std::shared_ptr<DOFHandlerMultiDim> dh_par = std::make_shared<DOFHandlerMultiDim>( const_cast<Mesh &>(*mesh) );
    std::shared_ptr<DiscreteSpace> ds = std::make_shared<EqualOrderDiscreteSpace>( &const_cast<Mesh &>(*mesh), fe_);
	dh_par->distribute_dofs(ds);
	dh_ = dh_par;
    unsigned int ndofs = dh_->max_elem_dofs();

	if (this->boundary_domain_) fill_boundary_dofs(); // temporary solution for boundary mesh
	else data_vec_ = VectorMPI::sequential( dh_->lsize() ); // allocate data_vec_

	// initialization data of value handlers
	FEValueInitData init_data;
	init_data.dh = dh_;
	init_data.data_vec = data_vec_;
	init_data.ndofs = ndofs;
	init_data.n_comp = this->n_comp();
	init_data.comp_index = 0;

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
		bool boundary;
		if (is_native || this->interpolation_==DataInterpolation::identic_msh || this->interpolation_==DataInterpolation::equivalent_msh) {
			n_entities = dh_->mesh()->n_elements();
			boundary = this->boundary_domain_;
		} else {
			n_entities = ReaderCache::get_mesh(reader_file_)->n_elements();
			boundary = false;
		}
		auto input_data_cache = ReaderCache::get_reader(reader_file_)->template get_element_data<double>(n_entities, n_components,
				boundary, this->component_idx_);
		CheckResult checked_data = ReaderCache::get_reader(reader_file_)->scale_and_check_limits(field_name_,
				this->unit_conversion_coefficient_, default_value_);


	    if (checked_data == CheckResult::not_a_number) {
	        THROW( ExcUndefElementValue() << EI_Field(field_name_) );
	    }

		if (is_native) {
			this->calculate_native_values(input_data_cache);
		} else if (this->interpolation_==DataInterpolation::identic_msh || this->interpolation_==DataInterpolation::equivalent_msh) {
			this->calculate_elementwise_values(input_data_cache);
		} else if (this->interpolation_==DataInterpolation::gauss_p0) {
			this->interpolate_gauss(input_data_cache);
		} else { // DataInterpolation::interp_p0
			this->interpolate_intersection(input_data_cache);
		}

		return true;
	} else return false;

}


template <int spacedim, class Value>
void FieldFE<spacedim, Value>::interpolate_gauss(ElementDataCache<double>::ComponentDataPtr data_vec)
{
	static const unsigned int quadrature_order = 4; // parameter of quadrature
	std::shared_ptr<Mesh> source_mesh = ReaderCache::get_mesh(reader_file_);
	std::vector<unsigned int> searched_elements; // stored suspect elements in calculating the intersection
	std::vector<arma::vec::fixed<3>> q_points; // real coordinates of quadrature points
	std::vector<double> q_weights; // weights of quadrature points
	unsigned int quadrature_size=0; // size of quadrature point and weight vector
	std::vector<double> sum_val(dh_->max_elem_dofs()); // sum of value of one quadrature point
	unsigned int elem_count; // count of intersect (source) elements of one quadrature point
	std::vector<double> elem_value(dh_->max_elem_dofs()); // computed value of one (target) element
	bool contains; // sign if source element contains quadrature point

	{
		// set size of vectors to maximal count of quadrature points
		QGauss quad(3, quadrature_order);
		q_points.resize(quad.size());
		q_weights.resize(quad.size());
	}

	Mesh *mesh;
	if (this->boundary_domain_) mesh = dh_->mesh()->get_bc_mesh();
	else mesh = dh_->mesh();
	for (auto cell : dh_->own_range()) {
		auto ele = cell.elm();
		std::fill(elem_value.begin(), elem_value.end(), 0.0);
		switch (cell.dim()) {
		case 0:
			quadrature_size = 1;
			q_points[0] = *ele.node(0);
			q_weights[0] = 1.0;
			break;
		case 1:
			quadrature_size = value_handler1_.compute_quadrature(q_points, q_weights, ele, quadrature_order);
			break;
		case 2:
			quadrature_size = value_handler2_.compute_quadrature(q_points, q_weights, ele, quadrature_order);
			break;
		case 3:
			quadrature_size = value_handler3_.compute_quadrature(q_points, q_weights, ele, quadrature_order);
			break;
		}
		searched_elements.clear();
		source_mesh->get_bih_tree().find_bounding_box(ele.bounding_box(), searched_elements);

		for (unsigned int i=0; i<quadrature_size; ++i) {
			std::fill(sum_val.begin(), sum_val.end(), 0.0);
			elem_count = 0;
			for (std::vector<unsigned int>::iterator it = searched_elements.begin(); it!=searched_elements.end(); it++) {
				ElementAccessor<3> elm = source_mesh->element_accessor(*it);
				contains=false;
				switch (elm->dim()) {
				case 0:
					contains = arma::norm(*elm.node(0) - q_points[i], 2) < 4*std::numeric_limits<double>::epsilon();
					break;
				case 1:
					contains = MappingP1<1,3>::contains_point(q_points[i], elm);
					break;
				case 2:
					contains = MappingP1<2,3>::contains_point(q_points[i], elm);
					break;
				case 3:
					contains = MappingP1<3,3>::contains_point(q_points[i], elm);
					break;
				default:
					ASSERT(false).error("Invalid element dimension!");
				}
				if ( contains ) {
					// projection point in element
					unsigned int index = sum_val.size() * (*it);
					for (unsigned int j=0; j < sum_val.size(); j++) {
						sum_val[j] += (*data_vec)[index+j];
					}
					++elem_count;
				}
			}

			if (elem_count > 0) {
				for (unsigned int j=0; j < sum_val.size(); j++) {
					elem_value[j] += (sum_val[j] / elem_count) * q_weights[i];
				}
			}
		}

		LocDofVec loc_dofs;
		if (this->boundary_domain_) loc_dofs = value_handler1_.get_loc_dof_indices(cell.elm_idx());
		else loc_dofs = cell.get_loc_dof_indices();

		ASSERT_LE_DBG(loc_dofs.n_elem, elem_value.size());
		for (unsigned int i=0; i < elem_value.size(); i++) {
			ASSERT_LT_DBG( loc_dofs[i], (int)data_vec_.size());
			data_vec_[loc_dofs[i]] = elem_value[i] * this->unit_conversion_coefficient_;
		}
	}
}


template <int spacedim, class Value>
void FieldFE<spacedim, Value>::interpolate_intersection(ElementDataCache<double>::ComponentDataPtr data_vec)
{
	std::shared_ptr<Mesh> source_mesh = ReaderCache::get_mesh(reader_file_);
	std::vector<unsigned int> searched_elements; // stored suspect elements in calculating the intersection
	std::vector<double> value(dh_->max_elem_dofs());
	double total_measure;
	double measure = 0;

	Mesh *mesh;
	if (this->boundary_domain_) mesh = dh_->mesh()->get_bc_mesh();
	else mesh = dh_->mesh();
	for (auto elm : mesh->elements_range()) {
		if (elm.dim() == 3) {
			xprintf(Err, "Dimension of element in target mesh must be 0, 1 or 2! elm.idx() = %d\n", elm.idx());
		}

		double epsilon = 4* numeric_limits<double>::epsilon() * elm.measure();

		// gets suspect elements
		if (elm.dim() == 0) {
			searched_elements.clear();
			source_mesh->get_bih_tree().find_point(*elm.node(0), searched_elements);
		} else {
			BoundingBox bb = elm.bounding_box();
			searched_elements.clear();
			source_mesh->get_bih_tree().find_bounding_box(bb, searched_elements);
		}

		// set zero values of value object
		std::fill(value.begin(), value.end(), 0.0);
		total_measure=0.0;

		START_TIMER("compute_pressure");
		ADD_CALLS(searched_elements.size());


        for (std::vector<unsigned int>::iterator it = searched_elements.begin(); it!=searched_elements.end(); it++)
        {
            ElementAccessor<3> ele = source_mesh->element_accessor(*it);
            if (ele->dim() == 3) {
                // get intersection (set measure = 0 if intersection doesn't exist)
                switch (elm.dim()) {
                    case 0: {
                        arma::vec::fixed<3> real_point = *elm.node(0);
                        arma::mat::fixed<3, 4> elm_map = MappingP1<3,3>::element_map(ele);
                        arma::vec::fixed<4> unit_point = MappingP1<3,3>::project_real_to_unit(real_point, elm_map);

                        measure = (std::fabs(arma::sum( unit_point )-1) <= 1e-14
                                        && arma::min( unit_point ) >= 0)
                                            ? 1.0 : 0.0;
                        break;
                    }
                    case 1: {
                        IntersectionAux<1,3> is;
                        ComputeIntersection<1,3> CI(elm, ele, source_mesh.get());
                        CI.init();
                        CI.compute(is);

                        IntersectionLocal<1,3> ilc(is);
                        measure = ilc.compute_measure() * elm.measure();
                        break;
                    }
                    case 2: {
                        IntersectionAux<2,3> is;
                        ComputeIntersection<2,3> CI(elm, ele, source_mesh.get());
                        CI.init();
                        CI.compute(is);

                        IntersectionLocal<2,3> ilc(is);
                        measure = 2 * ilc.compute_measure() * elm.measure();
                        break;
                    }
                }

				//adds values to value_ object if intersection exists
				if (measure > epsilon) {
					unsigned int index = value.size() * (*it);
			        std::vector<double> &vec = *( data_vec.get() );
			        for (unsigned int i=0; i < value.size(); i++) {
			        	value[i] += vec[index+i] * measure;
			        }
					total_measure += measure;
				}
			}
		}

		// computes weighted average, store it to data vector
		if (total_measure > epsilon) {
			VectorMPI::VectorDataPtr data_vector = data_vec_.data_ptr();

			LocDofVec loc_dofs;
			if (this->boundary_domain_) loc_dofs = value_handler1_.get_loc_dof_indices(elm.idx());
			else{
				DHCellAccessor cell = dh_->cell_accessor_from_element(elm.idx());
				loc_dofs = cell.get_loc_dof_indices();
			}

			ASSERT_LE_DBG(loc_dofs.n_elem, value.size());
			for (unsigned int i=0; i < value.size(); i++) {
				(*data_vector)[ loc_dofs[i] ] = value[i] / total_measure;
			}
		} else {
			WarningOut().fmt("Processed element with idx {} is out of source mesh!\n", elm.idx());
		}
		END_TIMER("compute_pressure");

	}
}


template <int spacedim, class Value>
void FieldFE<spacedim, Value>::calculate_native_values(ElementDataCache<double>::ComponentDataPtr data_cache)
{
	// Same algorithm as in output of Node_data. Possibly code reuse.
	unsigned int dof_size, data_vec_i;
	std::vector<unsigned int> count_vector(data_vec_.size(), 0);
	data_vec_.zero_entries();
	VectorMPI::VectorDataPtr data_vector = data_vec_.data_ptr();
	std::vector<LongIdx> global_dof_indices(dh_->max_elem_dofs());

	// iterate through cells, assembly MPIVector
	for (auto cell : dh_->own_range()) {
		dof_size = cell.get_dof_indices(global_dof_indices);
		LocDofVec loc_dofs = cell.get_loc_dof_indices();
		data_vec_i = cell.elm_idx() * dof_size;
		ASSERT_EQ_DBG(dof_size, loc_dofs.n_elem);
		for (unsigned int i=0; i<dof_size; ++i, ++data_vec_i) {
			(*data_vector)[ loc_dofs[i] ] += (*data_cache)[ global_dof_indices[i] ];
			++count_vector[ loc_dofs[i] ];
		}
	}

	// compute averages of values
	for (unsigned int i=0; i<data_vec_.size(); ++i) {
		if (count_vector[i]>0) (*data_vector)[i] /= count_vector[i];
	}
}


template <int spacedim, class Value>
void FieldFE<spacedim, Value>::calculate_elementwise_values(ElementDataCache<double>::ComponentDataPtr data_cache)
{
	// Same algorithm as in output of Node_data. Possibly code reuse.
	unsigned int data_vec_i;
	std::vector<unsigned int> count_vector(data_vec_.size(), 0);
	data_vec_.zero_entries();

	// iterate through elements, assembly global vector and count number of writes
	if (this->boundary_domain_) {
		Mesh *mesh = dh_->mesh()->get_bc_mesh();
		for (auto ele : mesh->elements_range()) { // remove special case for rank == 0 - necessary for correct output
			LocDofVec loc_dofs = value_handler1_.get_loc_dof_indices(ele.idx());
			data_vec_i = ele.idx() * dh_->max_elem_dofs();
			for (unsigned int i=0; i<loc_dofs.n_elem; ++i, ++data_vec_i) {
				ASSERT_LT_DBG(loc_dofs[i], data_vec_.size());
				data_vec_[ loc_dofs[i] ] += (*data_cache)[data_vec_i];
				++count_vector[ loc_dofs[i] ];
			}
		}
	}
	else {
		// iterate through cells, assembly global vector and count number of writes - prepared solution for further development
		for (auto cell : dh_->own_range()) {
			LocDofVec loc_dofs = cell.get_loc_dof_indices();
			data_vec_i = cell.elm_idx() * dh_->max_elem_dofs();
			for (unsigned int i=0; i<loc_dofs.n_elem; ++i, ++data_vec_i) {
				ASSERT_LT_DBG(loc_dofs[i], data_vec_.size());
				data_vec_[ loc_dofs[i] ] += (*data_cache)[data_vec_i];
				++count_vector[ loc_dofs[i] ];
			}
		}
	}

	// compute averages of values
	for (unsigned int i=0; i<data_vec_.size(); ++i) {
		if (count_vector[i]>0) data_vec_[i] /= count_vector[i];
	}
}


template <int spacedim, class Value>
void FieldFE<spacedim, Value>::native_data_to_cache(ElementDataCache<double> &output_data_cache) {
	ASSERT_EQ(output_data_cache.n_values() * output_data_cache.n_comp(), dh_->distr()->lsize()).error();
	double loc_values[output_data_cache.n_comp()];
	unsigned int i, dof_filled_size;

	VectorMPI::VectorDataPtr data_vec = data_vec_.data_ptr();
	for (auto dh_cell : dh_->own_range()) {
		LocDofVec loc_dofs = dh_cell.get_loc_dof_indices();
		for (i=0; i<loc_dofs.n_elem; ++i) loc_values[i] = (*data_vec)[ loc_dofs[i] ];
		for ( ; i<output_data_cache.n_comp(); ++i) loc_values[i] = numeric_limits<double>::signaling_NaN();
		output_data_cache.store_value( dh_cell.local_idx(), loc_values );
	}

	output_data_cache.set_dof_handler_hash( dh_->hash() );
}



template <int spacedim, class Value>
inline unsigned int FieldFE<spacedim, Value>::data_size() const {
	return data_vec_.size();
}



template <int spacedim, class Value>
void FieldFE<spacedim, Value>::local_to_ghost_data_scatter_begin() {
	data_vec_.local_to_ghost_begin();
}



template <int spacedim, class Value>
void FieldFE<spacedim, Value>::local_to_ghost_data_scatter_end() {
	data_vec_.local_to_ghost_end();
}



template <int spacedim, class Value>
FieldFE<spacedim, Value>::~FieldFE()
{}


// Instantiations of FieldFE
INSTANCE_ALL(FieldFE)
