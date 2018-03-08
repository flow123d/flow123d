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
 * @file    field_interpolated_p0.cc
 * @brief   
 */


#include "field_interpolated_p0.hh"
#include "fields/field_instances.hh"	// for instantiation macros
#include "system/system.hh"
#include "io/msh_gmshreader.h"
#include "mesh/bih_tree.hh"
#include "io/reader_cache.hh"
#include "system/sys_profiler.hh"

#include "fem/mapping_p1.hh"

#include "intersection/intersection_aux.hh"
#include "intersection/intersection_local.hh"
#include "intersection/compute_intersection.hh"

namespace it = Input::Type;

FLOW123D_FORCE_LINK_IN_CHILD(field_interpolated)



template <int spacedim, class Value>
const Input::Type::Record & FieldInterpolatedP0<spacedim, Value>::get_input_type()
{
    return it::Record("FieldInterpolatedP0", FieldAlgorithmBase<spacedim,Value>::template_name()+" Field interpolated from external mesh data and piecewise constant on mesh elements.")
        .derive_from(FieldAlgorithmBase<spacedim, Value>::get_input_type())
        .copy_keys(FieldAlgorithmBase<spacedim, Value>::get_field_algo_common_keys())
        .declare_key("mesh_data_file", IT::FileName::input(), IT::Default::obligatory(),
                "Input file with ASCII GMSH file format.")
        .declare_key("field_name", IT::String(), IT::Default::obligatory(),
                "The values of the Field are read from the ```$ElementData``` section with field name given by this key.")
		//.declare_key("unit", FieldAlgorithmBase<spacedim, Value>::get_input_type_unit_si(), it::Default::optional(),
		//		"Definition of unit.")
        .close();
}



template <int spacedim, class Value>
const int FieldInterpolatedP0<spacedim, Value>::registrar =
		Input::register_class< FieldInterpolatedP0<spacedim, Value>, unsigned int >("FieldInterpolatedP0") +
		FieldInterpolatedP0<spacedim, Value>::get_input_type().size();


template <int spacedim, class Value>
FieldInterpolatedP0<spacedim, Value>::FieldInterpolatedP0(const unsigned int n_comp)
: FieldAlgorithmBase<spacedim, Value>(n_comp)
{}



template <int spacedim, class Value>
void FieldInterpolatedP0<spacedim, Value>::init_from_input(const Input::Record &rec, const struct FieldAlgoBaseInitData& init_data) {
	this->init_unit_conversion_coefficient(rec, init_data);


	// read mesh, create tree
    {
       reader_file_ = FilePath( rec.val<FilePath>("mesh_data_file") );
       source_mesh_ = ReaderCache::get_mesh(reader_file_ );
	   // no call to mesh->setup_topology, we need only elements, no connectivity
    }
	bih_tree_ = new BIHTree( source_mesh_.get() );

    // allocate data_
	unsigned int data_size = source_mesh_->element.size() * (this->value_.n_rows() * this->value_.n_cols());
	data_ = std::make_shared<std::vector<typename Value::element_type>>();
	data_->resize(data_size);

	field_name_ = rec.val<std::string>("field_name");
}




template <int spacedim, class Value>
bool FieldInterpolatedP0<spacedim, Value>::set_time(const TimeStep &time) {
	OLD_ASSERT(source_mesh_, "Null mesh pointer of elementwise field: %s, did you call init_from_input(Input::Record)?\n", field_name_.c_str());
    if ( reader_file_ == FilePath() ) return false;
    
    //walkaround for the steady time governor - there is no data to be read in time==infinity
    //TODO: is it possible to check this before calling set_time?
    //if (time == numeric_limits< double >::infinity()) return false;
    
    // value of last computed element must be recalculated if time is changed
    computed_elm_idx_ = numeric_limits<unsigned int>::max();

    bool boundary_domain_ = false;
    BaseMeshReader::HeaderQuery header_query(field_name_, time.end(), OutputTime::DiscreteSpace::ELEM_DATA);
    ReaderCache::get_reader(reader_file_ )->find_header(header_query);
    data_ = ReaderCache::get_reader(reader_file_ )->template get_element_data<typename Value::element_type>(
    		source_mesh_->element.size(), this->value_.n_rows() * this->value_.n_cols(), boundary_domain_, this->component_idx_);
    this->scale_data();

    return true;
}



template <int spacedim, class Value>
typename Value::return_type const &FieldInterpolatedP0<spacedim, Value>::value(const Point &p, const ElementAccessor<spacedim> &elm)
{
	OLD_ASSERT( elm.is_elemental(), "FieldInterpolatedP0 works only for 'elemental' ElementAccessors.\n");
	if (elm.idx() != computed_elm_idx_ || elm.is_boundary() != computed_elm_boundary_) {
		computed_elm_idx_ = elm.idx();
		computed_elm_boundary_ = elm.is_boundary();

		if (elm.dim() == 3) {
			xprintf(Err, "Dimension of element in target mesh must be 0, 1 or 2! elm.idx() = %d\n", elm.idx());
		}

		double epsilon = 4* numeric_limits<double>::epsilon() * elm.element()->measure();

		// gets suspect elements
		if (elm.dim() == 0) {
			searched_elements_.clear();
			((BIHTree *)bih_tree_)->find_point(elm.element()->node[0]->point(), searched_elements_);
		} else {
			BoundingBox bb;
			elm.element()->get_bounding_box(bb);
			searched_elements_.clear();
			((BIHTree *)bih_tree_)->find_bounding_box(bb, searched_elements_);
		}

		// set zero values of value_ object
		for (unsigned int i=0; i < this->value_.n_rows(); i++) {
			for (unsigned int j=0; j < this->value_.n_cols(); j++) {
				this->value_(i,j) = 0.0;
			}
		}

		double total_measure=0.0, measure;
// 		ngh::TIntersectionType iType;

		START_TIMER("compute_pressure");
		ADD_CALLS(searched_elements_.size());
                
                ElementFullIter elm_full_iter = elm.full_iter();
                
                MappingP1<3,3> mapping;
                
		for (std::vector<unsigned int>::iterator it = searched_elements_.begin(); it!=searched_elements_.end(); it++)
		{
                    ElementFullIter ele = source_mesh_->element( *it );
                    if (ele->dim() == 3) {
                        // get intersection (set measure = 0 if intersection doesn't exist)
                        switch (elm.dim()) {
                            case 0: {
                                arma::vec::fixed<3> real_point = elm_full_iter->node[0]->point();
                                arma::mat::fixed<3, 4> elm_map = mapping.element_map(*elm_full_iter);
                                arma::vec::fixed<4> unit_point = mapping.project_real_to_unit(real_point, elm_map);
                                
                                measure = (std::fabs(arma::sum( unit_point )-1) <= 1e-14
                                                && arma::min( unit_point ) >= 0)
                                                    ? 1.0 : 0.0;
                                break;
                            }
                            case 1: {
                                IntersectionAux<1,3> is;
                                ComputeIntersection<1,3> CI(elm_full_iter, ele, source_mesh_.get());
                                CI.init();
                                CI.compute(is);

                                IntersectionLocal<1,3> ilc(is);
                                measure = ilc.compute_measure() * elm_full_iter->measure();
                                break;
                            }
                            case 2: {
                                IntersectionAux<2,3> is;
                                ComputeIntersection<2,3> CI(elm_full_iter, ele, source_mesh_.get());
                                CI.init();
                                CI.compute(is);

                                IntersectionLocal<2,3> ilc(is);
                                measure = 2 * ilc.compute_measure() * elm_full_iter->measure();
                                break;
                            }
                        }

                        //adds values to value_ object if intersection exists
                        if (measure > epsilon) {
                                unsigned int index = this->value_.n_rows() * this->value_.n_cols() * (*it);
                        std::vector<typename Value::element_type> &vec = *( data_.get() );
                        typename Value::return_type & ret_type_value = const_cast<typename Value::return_type &>( Value::from_raw(this->r_value_,  (typename Value::element_type *)(&vec[index])) );
                                Value tmp_value = Value( ret_type_value );

                                for (unsigned int i=0; i < this->value_.n_rows(); i++) {
                                        for (unsigned int j=0; j < this->value_.n_cols(); j++) {
                                                this->value_(i,j) += tmp_value(i,j) * measure;
                                        }
                                }
                                total_measure += measure;
                        }
                    }
		}
                
		// computes weighted average
		if (total_measure > epsilon) {
			for (unsigned int i=0; i < this->value_.n_rows(); i++) {
				for (unsigned int j=0; j < this->value_.n_cols(); j++) {
					this->value_(i,j) /= total_measure;
				}
			}
		} else {
			WarningOut().fmt("Processed element with idx {} is out of source mesh!\n", elm.idx());
		}
		END_TIMER("compute_pressure");

	}
    return this->r_value_;
}



template <int spacedim, class Value>
void FieldInterpolatedP0<spacedim, Value>::value_list(const std::vector< Point >  &point_list, const ElementAccessor<spacedim> &elm,
                       std::vector<typename Value::return_type>  &value_list)
{
	OLD_ASSERT( elm.is_elemental(), "FieldInterpolatedP0 works only for 'elemental' ElementAccessors.\n");
    FieldAlgorithmBase<spacedim, Value>::value_list(point_list, elm, value_list);
}



template <int spacedim, class Value>
void FieldInterpolatedP0<spacedim, Value>::scale_data()
{
	if (Value::is_scalable()) {
		std::vector<typename Value::element_type> &vec = *( data_.get() );
		for(unsigned int i=0; i<vec.size(); ++i)
			vec[i] *= this->unit_conversion_coefficient_;
	}
}



// Instantiations of FieldInterpolatedP0
INSTANCE_ALL(FieldInterpolatedP0)
