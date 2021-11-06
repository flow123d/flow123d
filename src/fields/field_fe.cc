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
#include "tools/unit_converter.hh"
#include "intersection/intersection_aux.hh"
#include "intersection/intersection_local.hh"
#include "intersection/compute_intersection.hh"




/// Implementation.

namespace it = Input::Type;




FLOW123D_FORCE_LINK_IN_CHILD(field_fe)



/************************************************************************************
 * Implementation of FieldFE methods
 */
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
        .declare_key("time_unit", UnitConverter::get_input_type(), TimeUnitConversion::get_input_default(),
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
  dh_(nullptr), field_name_(""), discretization_(OutputTime::DiscreteSpace::UNDEFINED), fe_values_(4)
{
	this->is_constant_in_space_ = false;
}


template<int spacedim, class Value>
typename Field<spacedim,Value>::FieldBasePtr FieldFE<spacedim, Value>::NativeFactory::create_field(Input::Record rec, const FieldCommon &field) {
	Input::Array multifield_arr;
	if (rec.opt_val(field.input_name(), multifield_arr))
	{
		unsigned int position = 0;
		auto it = multifield_arr.begin<Input::AbstractRecord>();
		if (multifield_arr.size() > 1)
			while (index_ != position) {
				++it; ++position;
			}

        Input::Record field_rec = (Input::Record)(*it);
        if (field_rec.val<std::string>("TYPE") == "FieldFE") {
            OutputTime::DiscreteSpace discretization;
            if (field_rec.opt_val<OutputTime::DiscreteSpace>("input_discretization", discretization)) {
                if (discretization == OutputTime::DiscreteSpace::NATIVE_DATA) {
                    std::shared_ptr< FieldFE<spacedim, Value> > field_fe = std::make_shared< FieldFE<spacedim, Value> >(field.n_comp());
                    FieldAlgoBaseInitData init_data(field.input_name(), field.n_comp(), field.units(), field.limits(), field.get_flags());
                    field_fe->init_from_input( *it, init_data );
                    field_fe->set_fe_data(conc_dof_handler_, dof_vector_);
                    return field_fe;
                }
            }
        }
	}

	return NULL;
}


template <int spacedim, class Value>
VectorMPI FieldFE<spacedim, Value>::set_fe_data(std::shared_ptr<DOFHandlerMultiDim> dh, VectorMPI dof_values, unsigned int block_index)
{
    dh_ = dh;
    if (dof_values.size()==0) { //create data vector according to dof handler - Warning not tested yet
        data_vec_ = dh_->create_vector();
        data_vec_.zero_entries();
    } else {
        data_vec_ = dof_values;
    }

    if ( block_index == FieldFE<spacedim, Value>::undef_uint ) {
        this->fill_fe_item<0>();
        this->fill_fe_item<1>();
        this->fill_fe_item<2>();
        this->fill_fe_item<3>();
        this->fe_ = dh_->ds()->fe();
    } else {
        this->fill_fe_system_data<0>(block_index);
        this->fill_fe_system_data<1>(block_index);
        this->fill_fe_system_data<2>(block_index);
        this->fill_fe_system_data<3>(block_index);
        this->fe_ = MixedPtr<FiniteElement>(
                std::dynamic_pointer_cast<FESystem<0>>( dh_->ds()->fe()[Dim<0>{}] )->fe()[block_index],
                std::dynamic_pointer_cast<FESystem<1>>( dh_->ds()->fe()[Dim<1>{}] )->fe()[block_index],
                std::dynamic_pointer_cast<FESystem<2>>( dh_->ds()->fe()[Dim<2>{}] )->fe()[block_index],
                std::dynamic_pointer_cast<FESystem<3>>( dh_->ds()->fe()[Dim<3>{}] )->fe()[block_index]
                );
    }

    unsigned int ndofs = dh_->max_elem_dofs();

    // initialization data of value handlers
	FEValueInitData init_data;
	init_data.dh = dh_;
	init_data.data_vec = data_vec_;
	init_data.ndofs = ndofs;
	init_data.n_comp = this->n_comp();
	init_data.mixed_fe = this->fe_;

	// initialize value handler objects
	init_data.range_begin = this->fe_item_[0].range_begin_;
	init_data.range_end = this->fe_item_[0].range_end_;
	value_handler0_.initialize(init_data);
	init_data.range_begin = this->fe_item_[1].range_begin_;
	init_data.range_end = this->fe_item_[1].range_end_;
	value_handler1_.initialize(init_data);
	init_data.range_begin = this->fe_item_[2].range_begin_;
	init_data.range_end = this->fe_item_[2].range_end_;
	value_handler2_.initialize(init_data);
	init_data.range_begin = this->fe_item_[3].range_begin_;
	init_data.range_end = this->fe_item_[3].range_end_;
	value_handler3_.initialize(init_data);

	// set interpolation
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
void FieldFE<spacedim, Value>::cache_update(FieldValueCache<typename Value::element_type> &data_cache,
		ElementCacheMap &cache_map, unsigned int region_patch_idx)
{
    ASSERT( !boundary_dofs_ ).error("boundary field NOT supported!!\n");
    Armor::ArmaMat<typename Value::element_type, Value::NRows_, Value::NCols_> mat_value;

    unsigned int reg_chunk_begin = cache_map.region_chunk_begin(region_patch_idx);
    unsigned int reg_chunk_end = cache_map.region_chunk_end(region_patch_idx);
    unsigned int last_element_idx = -1;
    DHCellAccessor cell = *( dh_->local_range().begin() ); //needs set variable for correct compiling
    LocDofVec loc_dofs;
    unsigned int range_bgn=0, range_end=0;

    for (unsigned int i_data = reg_chunk_begin; i_data < reg_chunk_end; ++i_data) { // i_eval_point_data
        unsigned int elm_idx = cache_map.eval_point_data(i_data).i_element_;
        if (elm_idx != last_element_idx) {
            ElementAccessor<spacedim> elm(dh_->mesh(), elm_idx);
            fe_values_[elm.dim()].reinit( elm );
            cell = dh_->cell_accessor_from_element( elm_idx );
            loc_dofs = cell.get_loc_dof_indices();
            last_element_idx = elm_idx;
            range_bgn = this->fe_item_[elm.dim()].range_begin_;
            range_end = this->fe_item_[elm.dim()].range_end_;
        }

        unsigned int i_ep=cache_map.eval_point_data(i_data).i_eval_point_;
        //DHCellAccessor cache_cell = cache_map(cell);
        mat_value.fill(0.0);
        for (unsigned int i_dof=range_bgn, i_cdof=0; i_dof<range_end; i_dof++, i_cdof++) {
            mat_value += data_vec_.get(loc_dofs[i_dof]) * this->handle_fe_shape(cell.dim(), i_cdof, i_ep);
        }
        data_cache.set(i_data) = mat_value;
    }
}


template <int spacedim, class Value>
void FieldFE<spacedim, Value>::cache_reinit(const ElementCacheMap &cache_map)
{
    std::shared_ptr<EvalPoints> eval_points = cache_map.eval_points();
    std::array<Quadrature, 4> quads{QGauss(0, 1), this->init_quad<1>(eval_points), this->init_quad<2>(eval_points), this->init_quad<3>(eval_points)};
    fe_values_[0].initialize(quads[0], *this->fe_[0_d], update_values);
    fe_values_[1].initialize(quads[1], *this->fe_[1_d], update_values);
    fe_values_[2].initialize(quads[2], *this->fe_[2_d], update_values);
    fe_values_[3].initialize(quads[3], *this->fe_[3_d], update_values);
}


template <int spacedim, class Value>
template <unsigned int dim>
Quadrature FieldFE<spacedim, Value>::init_quad(std::shared_ptr<EvalPoints> eval_points)
{
    Quadrature quad(dim, eval_points->size(dim));
    for (unsigned int k=0; k<eval_points->size(dim); k++)
        quad.set(k) = eval_points->local_point<dim>(k);
    return quad;
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
        if (this->interpolation_ == DataInterpolation::identic_msh) {
            ReaderCache::get_element_ids(reader_file_, *mesh);
        } else {
            auto source_mesh = ReaderCache::get_mesh(reader_file_);
            ReaderCache::get_element_ids(reader_file_, *(source_mesh.get()));
            if (this->interpolation_ == DataInterpolation::equivalent_msh) {
                source_target_mesh_elm_map_ = ReaderCache::get_target_mesh_element_map(reader_file_, const_cast<Mesh *>(mesh));
                if (source_target_mesh_elm_map_->empty()) { // incompatible meshes
                    this->interpolation_ = DataInterpolation::gauss_p0;
                    WarningOut().fmt("Source mesh of FieldFE '{}' is not compatible with target mesh.\nInterpolation of input data will be changed to 'P0_gauss'.\n",
                            field_name_);
                }
            } else if (this->interpolation_ == DataInterpolation::interp_p0) {
                if (!boundary_domain) {
                    this->interpolation_ = DataInterpolation::gauss_p0;
                    WarningOut().fmt("Interpolation 'P0_intersection' of FieldFE '{}' can't be used on bulk region.\nIt will be changed to 'P0_gauss'.\n",
                            field_name_);
                }
            }
        }
        if (dh_ == nullptr) this->make_dof_handler(mesh);
	}
}



template <int spacedim, class Value>
void FieldFE<spacedim, Value>::fill_boundary_dofs() {
	ASSERT(this->boundary_domain_);

	auto bc_mesh = dh_->mesh()->get_bc_mesh();
	unsigned int n_comp = this->value_.n_rows() * this->value_.n_cols();
	boundary_dofs_ = std::make_shared< std::vector<IntIdx> >( n_comp * bc_mesh->n_elements() );
	std::vector<IntIdx> &in_vec = *( boundary_dofs_.get() );
	unsigned int j = 0; // actual index to boundary_dofs_ vector

	for (auto ele : bc_mesh->elements_range()) {
		IntIdx elm_shift = n_comp * ele.idx();
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
	MixedPtr<FiniteElement> fe;
	switch (this->value_.n_rows() * this->value_.n_cols()) { // by number of components
		case 1: { // scalar
			fe = MixedPtr<FE_P_disc>(0);
			break;
		}
		case 3: { // vector
			 MixedPtr<FE_P_disc>   fe_base(0) ;
			fe = mixed_fe_system(fe_base, FEType::FEVector, 3);
			break;
		}
		case 9: { // tensor
		    MixedPtr<FE_P_disc>   fe_base(0) ;
            fe = mixed_fe_system(fe_base, FEType::FETensor, 9);
			break;
		}
		default:
			ASSERT(false).error("Should not happen!\n");
	}

	std::shared_ptr<DOFHandlerMultiDim> dh_par = std::make_shared<DOFHandlerMultiDim>( const_cast<Mesh &>(*mesh) );
    std::shared_ptr<DiscreteSpace> ds = std::make_shared<EqualOrderDiscreteSpace>( &const_cast<Mesh &>(*mesh), fe);
	dh_par->distribute_dofs(ds);
	dh_ = dh_par;
    unsigned int ndofs = dh_->max_elem_dofs();

    this->fill_fe_item<0>();
    this->fill_fe_item<1>();
    this->fill_fe_item<2>();
    this->fill_fe_item<3>();
    this->fe_ = dh_->ds()->fe();

    if (this->boundary_domain_) fill_boundary_dofs(); // temporary solution for boundary mesh
    else data_vec_ = VectorMPI::sequential( dh_->lsize() ); // allocate data_vec_

	// initialization data of value handlers
	FEValueInitData init_data;
	init_data.dh = dh_;
	init_data.data_vec = data_vec_;
	init_data.ndofs = ndofs;
	init_data.n_comp = this->n_comp();
	init_data.mixed_fe = this->fe_;

	// initialize value handler objects
	init_data.range_begin = this->fe_item_[0].range_begin_;
	init_data.range_end = this->fe_item_[0].range_end_;
	value_handler0_.initialize(init_data);
	init_data.range_begin = this->fe_item_[1].range_begin_;
	init_data.range_end = this->fe_item_[1].range_end_;
	value_handler1_.initialize(init_data);
	init_data.range_begin = this->fe_item_[2].range_begin_;
	init_data.range_end = this->fe_item_[2].range_end_;
	value_handler2_.initialize(init_data);
	init_data.range_begin = this->fe_item_[3].range_begin_;
	init_data.range_end = this->fe_item_[3].range_end_;
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
		double time_unit_coef = time.read_coef(in_rec_.find<Input::Record>("time_unit"));
		double time_shift = time.read_time( in_rec_.find<Input::Tuple>("read_time_shift") );
		double read_time = (time.end()+time_shift) / time_unit_coef;
		BaseMeshReader::HeaderQuery header_query(field_name_, read_time, this->discretization_, dh_->hash());
		ReaderCache::get_reader(reader_file_)->find_header(header_query);
		// TODO: use default and check NaN values in data_vec

		unsigned int n_entities;
		bool is_native = (header_query.discretization == OutputTime::DiscreteSpace::NATIVE_DATA);
		bool boundary;
		if (is_native || this->interpolation_==DataInterpolation::identic_msh || this->interpolation_==DataInterpolation::equivalent_msh) {
			boundary = this->boundary_domain_;
		} else {
			boundary = false;
		}
		if (is_native) {
		    n_entities = boundary ? dh_->mesh()->get_bc_mesh()->n_elements() : dh_->mesh()->n_elements();
		    n_components *= dh_->max_elem_dofs();
		} else if (this->interpolation_==DataInterpolation::identic_msh) {
			n_entities = boundary ? dh_->mesh()->get_bc_mesh()->n_elements() : dh_->mesh()->n_elements();
		} else {
			n_entities = boundary ? ReaderCache::get_mesh(reader_file_)->get_bc_mesh()->n_elements() : ReaderCache::get_mesh(reader_file_)->n_elements();
		}
		auto input_data_cache = ReaderCache::get_reader(reader_file_)->template get_element_data<double>(n_entities, n_components,
				boundary, this->component_idx_);
		CheckResult checked_data = ReaderCache::get_reader(reader_file_)->scale_and_check_limits(field_name_,
				this->unit_conversion_coefficient_, default_value_);


	    if ( !is_native && (checked_data == CheckResult::not_a_number) ) {
	        THROW( ExcUndefElementValue() << EI_Field(field_name_) );
	    }

		if (is_native) {
			this->calculate_native_values(input_data_cache);
		} else if (this->interpolation_==DataInterpolation::identic_msh) {
			this->calculate_identic_values(input_data_cache);
		} else if (this->interpolation_==DataInterpolation::equivalent_msh) {
			this->calculate_equivalent_values(input_data_cache);
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
			data_vec_.set( loc_dofs[i], elem_value[i] * this->unit_conversion_coefficient_ );
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
			THROW( ExcInvalidElemeDim() << EI_ElemIdx(elm.idx()) );
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
            ElementAccessor<3> source_elm = source_mesh->element_accessor(*it);
            if (source_elm->dim() == 3) {
                // get intersection (set measure = 0 if intersection doesn't exist)
                switch (elm.dim()) {
                    case 0: {
                        arma::vec::fixed<3> real_point = *elm.node(0);
                        arma::mat::fixed<3, 4> elm_map = MappingP1<3,3>::element_map(source_elm);
                        arma::vec::fixed<4> unit_point = MappingP1<3,3>::project_real_to_unit(real_point, elm_map);

                        measure = (std::fabs(arma::sum( unit_point )-1) <= 1e-14
                                        && arma::min( unit_point ) >= 0)
                                            ? 1.0 : 0.0;
                        break;
                    }
                    case 1: {
                        IntersectionAux<1,3> is(elm.mesh_idx(), source_elm.mesh_idx());
                        ComputeIntersection<1,3> CI(elm, source_elm, source_mesh.get());
                        CI.init();
                        CI.compute(is);

                        IntersectionLocal<1,3> ilc(is);
                        measure = ilc.compute_measure() * elm.measure();
                        break;
                    }
                    case 2: {
                        IntersectionAux<2,3> is(elm.mesh_idx(), source_elm.mesh_idx());
                        ComputeIntersection<2,3> CI(elm, source_elm, source_mesh.get());
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
			LocDofVec loc_dofs;
			if (this->boundary_domain_) loc_dofs = value_handler1_.get_loc_dof_indices(elm.idx());
			else{
				DHCellAccessor cell = dh_->cell_accessor_from_element(elm.idx());
				loc_dofs = cell.get_loc_dof_indices();
			}

			ASSERT_LE_DBG(loc_dofs.n_elem, value.size());
			for (unsigned int i=0; i < value.size(); i++) {
				data_vec_.set(loc_dofs[i], value[i] / total_measure);
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
	std::vector<LongIdx> global_dof_indices(dh_->max_elem_dofs());
	std::vector<LongIdx> &source_target_vec = source_target_mesh_elm_map_->bulk;

	// iterate through cells, assembly MPIVector
	for (auto cell : dh_->own_range()) {
		dof_size = cell.get_dof_indices(global_dof_indices);
		LocDofVec loc_dofs = cell.get_loc_dof_indices();
		data_vec_i = source_target_vec[cell.elm_idx()] * dof_size;
		ASSERT_EQ_DBG(dof_size, loc_dofs.n_elem);
		for (unsigned int i=0; i<dof_size; ++i, ++data_vec_i) {
		    data_vec_.add( loc_dofs[i], (*data_cache)[ data_vec_i ] );
		    ++count_vector[ loc_dofs[i] ];
		}
	}

	// compute averages of values
	for (unsigned int i=0; i<data_vec_.size(); ++i) {
		if (count_vector[i]>0) data_vec_.normalize(i, count_vector[i]);
	}
}


template <int spacedim, class Value>
void FieldFE<spacedim, Value>::calculate_identic_values(ElementDataCache<double>::ComponentDataPtr data_cache)
{
	// Same algorithm as in output of Node_data. Possibly code reuse.
	unsigned int data_vec_i, i_elm;
	std::vector<unsigned int> count_vector(data_vec_.size(), 0);
	data_vec_.zero_entries();

	if (this->boundary_domain_) {
		// iterate through elements, assembly global vector and count number of writes
		Mesh *mesh = dh_->mesh()->get_bc_mesh();
		i_elm=0;
		for (auto ele : mesh->elements_range()) {
			LocDofVec loc_dofs = value_handler1_.get_loc_dof_indices(ele.idx());
			data_vec_i = i_elm * dh_->max_elem_dofs();
			for (unsigned int i=0; i<loc_dofs.n_elem; ++i, ++data_vec_i) {
				ASSERT_LT_DBG(loc_dofs[i], (LongIdx)data_vec_.size());
				data_vec_.add( loc_dofs[i], (*data_cache)[data_vec_i] );
				++count_vector[ loc_dofs[i] ];
			}
			i_elm++;
		}
	}
	else {
		// iterate through cells, assembly global vector and count number of writes - prepared solution for further development
		i_elm=0;
		for (auto cell : dh_->own_range()) {
			LocDofVec loc_dofs = cell.get_loc_dof_indices();
			data_vec_i = i_elm * dh_->max_elem_dofs();
			for (unsigned int i=0; i<loc_dofs.n_elem; ++i, ++data_vec_i) {
				ASSERT_LT_DBG(loc_dofs[i], (LongIdx)data_vec_.size());
				data_vec_.add( loc_dofs[i], (*data_cache)[data_vec_i] );
				++count_vector[ loc_dofs[i] ];
			}
			i_elm++;
		}
	}

	// compute averages of values
	for (unsigned int i=0; i<data_vec_.size(); ++i) {
		if (count_vector[i]>0) data_vec_.normalize(i, count_vector[i]);
	}
}


template <int spacedim, class Value>
void FieldFE<spacedim, Value>::calculate_equivalent_values(ElementDataCache<double>::ComponentDataPtr data_cache)
{
	// Same algorithm as in output of Node_data. Possibly code reuse.
	unsigned int data_vec_i;
	std::vector<unsigned int> count_vector(data_vec_.size(), 0);
	data_vec_.zero_entries();

	// iterate through elements, assembly global vector and count number of writes
	if (this->boundary_domain_) {
        std::vector<LongIdx> &source_target_vec = source_target_mesh_elm_map_->boundary;
		Mesh *mesh = dh_->mesh()->get_bc_mesh();
		for (auto ele : mesh->elements_range()) {
			LocDofVec loc_dofs = value_handler1_.get_loc_dof_indices(ele.idx());
			if (source_target_vec[ele.idx()] == (int)(undef_idx)) { // undefined value in input data mesh
				if ( std::isnan(default_value_) )
					THROW( ExcUndefElementValue() << EI_Field(field_name_) );
				for (unsigned int i=0; i<loc_dofs.n_elem; ++i) {
					ASSERT_LT_DBG(loc_dofs[i], (LongIdx)data_vec_.size());
					data_vec_.add( loc_dofs[i], default_value_ * this->unit_conversion_coefficient_ );
					++count_vector[ loc_dofs[i] ];
				}
			} else {
				data_vec_i = source_target_vec[ele.idx()] * dh_->max_elem_dofs();
				for (unsigned int i=0; i<loc_dofs.n_elem; ++i, ++data_vec_i) {
					ASSERT_LT_DBG(loc_dofs[i], (LongIdx)data_vec_.size());
					data_vec_.add( loc_dofs[i], (*data_cache)[data_vec_i] );
					++count_vector[ loc_dofs[i] ];
				}
			}
		}
	}
	else {
        std::vector<LongIdx> &source_target_vec = source_target_mesh_elm_map_->bulk;
		// iterate through cells, assembly global vector and count number of writes - prepared solution for further development
		for (auto cell : dh_->own_range()) {
			LocDofVec loc_dofs = cell.get_loc_dof_indices();
			if (source_target_vec[cell.elm_idx()] == (int)(undef_idx)) { // undefined value in input data mesh
				if ( std::isnan(default_value_) )
					THROW( ExcUndefElementValue() << EI_Field(field_name_) );
				for (unsigned int i=0; i<loc_dofs.n_elem; ++i) {
					ASSERT_LT_DBG(loc_dofs[i], (LongIdx)data_vec_.size());
					data_vec_.add( loc_dofs[i], default_value_ * this->unit_conversion_coefficient_ );
					++count_vector[ loc_dofs[i] ];
				}
			} else {
				data_vec_i = source_target_vec[cell.elm_idx()] * dh_->max_elem_dofs();
				for (unsigned int i=0; i<loc_dofs.n_elem; ++i, ++data_vec_i) {
					ASSERT_LT_DBG(loc_dofs[i], (LongIdx)data_vec_.size());
					data_vec_.add( loc_dofs[i], (*data_cache)[data_vec_i] );
					++count_vector[ loc_dofs[i] ];
				}
			}
		}
	}

	// compute averages of values
	for (unsigned int i=0; i<data_vec_.size(); ++i) {
		if (count_vector[i]>0) data_vec_.normalize(i, count_vector[i]);
	}
}


template <int spacedim, class Value>
void FieldFE<spacedim, Value>::native_data_to_cache(ElementDataCache<double> &output_data_cache) {
	//ASSERT_EQ(output_data_cache.n_values() * output_data_cache.n_comp(), dh_->distr()->lsize()).error();
	unsigned int n_vals = output_data_cache.n_comp() * output_data_cache.n_dofs_per_element();
	double loc_values[n_vals];
	unsigned int i;

	for (auto dh_cell : dh_->own_range()) {
		LocDofVec loc_dofs = dh_cell.get_loc_dof_indices();
		for (i=0; i<loc_dofs.n_elem; ++i) loc_values[i] = data_vec_.get( loc_dofs[i] );
		for ( ; i<n_vals; ++i) loc_values[i] = numeric_limits<double>::signaling_NaN();
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



/*template <int spacedim, class Value>
Armor::ArmaMat<typename Value::element_type, Value::NRows_, Value::NCols_> FieldFE<spacedim, Value>::handle_fe_shape(unsigned int dim,
        unsigned int i_dof, unsigned int i_qp, unsigned int comp_index)
{
    Armor::ArmaMat<typename Value::element_type, Value::NCols_, Value::NRows_> v;
    for (unsigned int c=0; c<Value::NRows_*Value::NCols_; ++c)
        v(c/spacedim,c%spacedim) = fe_values_[dim].shape_value_component(i_dof, i_qp, comp_index+c);
    if (Value::NRows_ == Value::NCols_)
        return v;
    else
        return v.t();
}*/



template <int spacedim, class Value>
FieldFE<spacedim, Value>::~FieldFE()
{}


// Instantiations of FieldFE
INSTANCE_ALL(FieldFE)
